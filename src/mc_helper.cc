#include "mc_helper.h"
#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include "constants.h"

namespace gagatt {
std::vector<ZCacheEntry> buildLumiCache(const MCConfig &cfg, double pc1,
                                        double pc2, double sqrts_min,
                                        double d_sqrts) {
    std::vector<ZCacheEntry> zcache(cfg.n_sqrts);

    const int print_every = std::max(1, cfg.n_sqrts / 5);
    std::cout << "-- precomputing lumi cache ...\n";
    for (int j = 0; j < cfg.n_sqrts; ++j) {
        const double sqrt_s_hat = sqrts_min + (j + 0.5) * d_sqrts;
        const double z = sqrt_s_hat / cfg.sqrt_s;
        const auto [lw, L_tot] =
            lumiWeightsAndTotal(z, cfg.x, cfg.pe1, pc1, cfg.pe2, pc2);
        zcache[j] = {lw, L_tot};

        if ((j + 1) % print_every == 0) {
            std::cout << std::format(" lumi cache: {}/{}\n", j + 1,
                                     cfg.n_sqrts);
        }
    }
    return zcache;
}

// -----------------------------------------------------------------------
// partialXsec:
// d sigma_hat / d cos_th (helicity-summed, luminosity-weighted)
// = (beta Nc / 32 pi s_hat) * |A_C|^2 * (C_1^w + C_3^w)
//
// sdc.norm_factor = C_1^w + C_3^w, with |A_C|^2 already absorbed via
// overall_fac^2 = (COUPLING_FACTOR / denom)^2 inside polCoeffsForHelicity.
// -----------------------------------------------------------------------
double partialXsec(double sqrt_s_hat, const SDMatrixCoefficients &sdc) {
    if (sdc.norm_factor <= 0.0) { return 0.0; }

    const double s_hat = sqrt_s_hat * sqrt_s_hat;
    const double r = 4.0 * MTOP2 / s_hat;
    if (r >= 1.0) { return 0.0; }

    const double beta = std::sqrt(1.0 - r);

    // beta Nc / (32 pi s_hat)
    const double prefac = beta * NC / (32.0 * std::numbers::pi * s_hat);
    return prefac * sdc.norm_factor;
}

// -----------------------------------------------------------------------
// eventRate:
// d^2 sigma / (d sqrt_s_hat d cos_th)
//
// = partialXsec
//   * L_tot(z) (photon lumi per unit z = sqrt_tau)
//   * GEV2_TO_FB
//   / sqrt_s (Jacobian: dz = d(sqrt_s_hat)/sqrt_s)
// -----------------------------------------------------------------------
double eventRate(double sqrt_s_hat, const SDMatrixCoefficients &sdc,
                 double L_tot, double sqrt_s) {
    const double xsec = partialXsec(sqrt_s_hat, sdc);
    if (xsec <= 0.0 || L_tot <= 0.0) { return 0.0; }

    return xsec * L_tot * GEV2_TO_FB / sqrt_s;
}

WeightTable buildWeightTable(const MCConfig &cfg,
                             const std::vector<ZCacheEntry> &zcache,
                             double sqrts_min, double d_sqrts, double d_cos) {
    WeightTable wt;
    const int N_bins = cfg.n_cos * cfg.n_sqrts;
    wt.bin_weights.assign(N_bins, 0.0);
    wt.sdc_cache.reserve(N_bins);

    double tw_con = 0.0, tw_D = 0.0, tw_m12 = 0.0;
    Eigen::Vector3d tw_bp = Eigen::Vector3d::Zero();
    Eigen::Vector3d tw_bm = Eigen::Vector3d::Zero();

    const int print_every = std::max(1, cfg.n_cos / 5);
    std::cout << "-- building weight table ...\n";
    for (int i = 0; i < cfg.n_cos; ++i) {
        const double cos_th = cfg.cos_th_min + (i + 0.5) * d_cos;
        for (int j = 0; j < cfg.n_sqrts; ++j) {
            const double sqrt_s_hat = sqrts_min + (j + 0.5) * d_sqrts;
            const SDMatrixCoefficients sdc(sqrt_s_hat, cos_th, zcache[j].lw);
            const Matrix4cd rho = spinDensityMatrix(sdc);

            const double rate =
                eventRate(sqrt_s_hat, sdc, zcache[j].L_tot, cfg.sqrt_s) *
                d_sqrts * d_cos;

            const int idx = i * cfg.n_sqrts + j;
            wt.bin_weights[idx] = std::max(0.0, rate);
            wt.sdc_cache.push_back(sdc);

            wt.total_weight += wt.bin_weights[idx];
            tw_con += wt.bin_weights[idx] * getConcurrence(rho);
            tw_D += wt.bin_weights[idx] * entanglementMarker(sdc);
            tw_m12 += wt.bin_weights[idx] * horodeckiMeasure(sdc);
            tw_bp += wt.bin_weights[idx] * sdc.bp;
            tw_bm += wt.bin_weights[idx] * sdc.bm;
        }
        if ((i + 1) % print_every == 0) {
            std::cout << std::format(" weight table: {}/{}\n", i + 1,
                                     cfg.n_cos);
        }
    }

    if (wt.total_weight <= 0.0) {
        std::cerr << "ERROR: total weight is zero — check kinematics.\n";
        wt.bin_weights.clear();
        return wt;
    }

    wt.theory_concurrence = tw_con / wt.total_weight;
    wt.theory_D = tw_D / wt.total_weight;
    wt.theory_m12 = tw_m12 / wt.total_weight;
    wt.theory_bp = tw_bp / wt.total_weight;
    wt.theory_bm = tw_bm / wt.total_weight;

    return wt;
}

// -----------------------------------------------------------------------
// sampleDecayAngles:
// Draw (q+, q-): unit vectors for l+ and l- in the top / anti-top
// rest frames from the joint angular distribution:
//
//   P(q+, q-) = (1/4pi)^2 * [1 + B+.q+ - B-.q- - q+.C.q-]
//
// Signs: +B+.q+  (top polarization)
//        -B-.q-  (anti-top polarization, minus sign)
//        -q+.C.q- (spin correlation, minus sign)
//
// Accept/reject with envelope w_max = 1 + |B+| + |B-| + ||C||_F.
// Returns the pair {q+, q-} as unit Eigen::Vector3d.
// -----------------------------------------------------------------------
std::pair<Eigen::Vector3d, Eigen::Vector3d> sampleDecayAngles(
    const SDMatrixCoefficients &sdc, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> uni_cos(-1.0, 1.0);
    std::uniform_real_distribution<double> uni_phi(0.0, 2.0 * std::numbers::pi);
    std::uniform_real_distribution<double> uni01(0.0, 1.0);

    // Conservative envelope (Frobenius norm upper-bounds |q+.C.q-|)
    const double w_max = 1.0 + sdc.bp.norm() + sdc.bm.norm() + sdc.cc.norm();

    for (;;) {
        // Sample q+ uniformly on S^2
        const double cth_p = uni_cos(rng);
        const double sth_p = std::sqrt(1.0 - cth_p * cth_p);
        const double phi_p = uni_phi(rng);
        const Eigen::Vector3d qp(sth_p * std::cos(phi_p),
                                 sth_p * std::sin(phi_p), cth_p);

        // Sample q- uniformly on S^2
        const double cth_m = uni_cos(rng);
        const double sth_m = std::sqrt(1.0 - cth_m * cth_m);
        const double phi_m = uni_phi(rng);
        const Eigen::Vector3d qm(sth_m * std::cos(phi_m),
                                 sth_m * std::sin(phi_m), cth_m);

        // Weight: 1 + B+.q+ - B-.q- - q+.C.q-
        const double w =
            1.0 + sdc.bp.dot(qp) - sdc.bm.dot(qm) - qp.dot(sdc.cc * qm);

        if (uni01(rng) * w_max <= w) { return {qp, qm}; }
    }
}

// -----------------------------------------------------------------------
// runEventLoop:
//
// For each event, draw a (cos_th, sqrt_s_hat) bin, then sample decay
// directions (q+, q-).
//
// Accumulate:
//   S_ij  = sum q+_i * q-_j          (first moment)
//   S2_ij = sum (q+_i * q-_j)^2      (second moment, for variance)
// -----------------------------------------------------------------------
EventLoopResult runEventLoop(const MCConfig &cfg, long long n_events,
                             const WeightTable &wt, std::mt19937_64 &rng,
                             bool verbose) {
    std::discrete_distribution<int> bin_dist(wt.bin_weights.begin(),
                                             wt.bin_weights.end());

    EventLoopResult ev;
    const long long print_every = std::max(1LL, n_events / 10);

    while (ev.n_accepted < n_events) {
        const int k = bin_dist(rng);
        const SDMatrixCoefficients &sdc = wt.sdc_cache[k];

        // Sample q+ and q- from the joint distribution
        const auto [qp, qm] = sampleDecayAngles(sdc, rng);

        // Accumulate first moments (for B_+, B_-)
        ev.S1_qp += qp;
        ev.S2_qp += qp.cwiseProduct(qp);
        ev.S1_qm += qm;
        ev.S2_qm += qm.cwiseProduct(qm);

        // Accumulate outer product and its element-wise square (for C_ij)
        const Eigen::Matrix3d outer = qp * qm.transpose();
        ev.S1_qpqm += outer;
        ev.S2_qpqm += outer.cwiseProduct(outer);

        ++ev.n_accepted;
        if (verbose && ev.n_accepted % print_every == 0) {
            std::cout << std::format(" events: {}/{}\n", ev.n_accepted,
                                     cfg.n_events);
        }
    }
    return ev;
}

// -----------------------------------------------------------------------
// reconstructFromMoments:
//
// After N events:
//   <q+_i>  = S1_qp/N,   B_{+i} = 3 <q+_i>
//   <q-_j>  = S1_qm/N,   B_{-j} = -3 <q-_j>
//   <qi qj> = S1_qpqm/N, C_ij   = -9 <q+_i q-_j>
//
//   Var[<q+_i>]   = (S2_qp/N - <q+_i>^2) / N,  sigma[B_+i] = 3 sqrt(Var)
//   Var[<q-_j>]   = (S2_qm/N - <q-_j>^2) / N,  sigma[B_-j] = 3 sqrt(Var)
//   Var[<qi qj>]  = (S2_qpqm/N - <qi qj>^2)/N, sigma[C_ij] = 9 sqrt(Var)
// -----------------------------------------------------------------------
ReconstructedMC reconstructFromMoments(const EventLoopResult &ev) {
    ReconstructedMC r;
    const double n = static_cast<double>(ev.n_accepted);

    // --- B_+ and B_- ---
    const Eigen::Vector3d mean_qp = ev.S1_qp / n;
    const Eigen::Vector3d mean_qp2 = ev.S2_qp / n;
    const Eigen::Vector3d var_qp =
        (mean_qp2 - mean_qp.cwiseProduct(mean_qp)) / n;

    r.mc_bp = 3.0 * mean_qp;
    r.sigma_bp = 3.0 * var_qp.cwiseMax(0.0).cwiseSqrt();

    const Eigen::Vector3d mean_qm = ev.S1_qm / n;
    const Eigen::Vector3d mean_qm2 = ev.S2_qm / n;
    const Eigen::Vector3d var_qm =
        (mean_qm2 - mean_qm.cwiseProduct(mean_qm)) / n;

    r.mc_bm = -3.0 * mean_qm;
    r.sigma_bm = 3.0 * var_qm.cwiseMax(0.0).cwiseSqrt();

    // --- C_ij ---
    const Eigen::Matrix3d mean_qpqm = ev.S1_qpqm / n;

    // C_ij^MC = -9 * <q+_i q-_j>
    r.mc_cij = -9.0 * mean_qpqm;

    // Variance of the mean per element: Var[<x>] = (E[x^2] - E[x]^2) / N
    const Eigen::Matrix3d mean_qpqm2 = ev.S2_qpqm / n;
    const Eigen::Matrix3d var_mean =
        (mean_qpqm2 - mean_qpqm.cwiseProduct(mean_qpqm)) / n;

    // sigma[C_ij] = 9 * sqrt(Var[<qi qj>])
    r.sigma_cij = 9.0 * var_mean.cwiseMax(0.0).cwiseSqrt();

    // Reconstruct density matrix from B+^MC, B-^MC, C_ij^MC
    const Matrix4cd mc_rho = reconstructRho(r.mc_bp, r.mc_bm, r.mc_cij);
    r.mc_concurrence = getConcurrence(mc_rho);
    // sigma_concurrence = (1/2) ||sigma_cij||_F
    // Uses all 9 elements of sigma_cij (Frobenius norm), which is conservative
    // for any C regardless of whether off-diagonal elements are active.
    // Exact propagation through the concurrence formula would require
    // eigenvectors of C*C^T; this bound is sufficient for significance
    // estimates.
    r.sigma_concurrence = 0.5 * r.sigma_cij.norm();
    r.significance_concurrence =
        (r.sigma_concurrence > 0.0 && r.mc_concurrence > 0.0)
            ? r.mc_concurrence / r.sigma_concurrence
            : 0.0;

    // D = (C_nn - |C_rr + C_kk|) / 3
    r.mc_D = entanglementMarker(r.mc_cij);
    // |dD/dC_nn| = |dD/dC_rr| = |dD/dC_kk| = 1/3 regardless of sign of
    // (C_rr+C_kk), so sigma_D has the same form as sigma_Tr[C]/3.
    const double sigma_tr_c =
        9.0 * std::sqrt(std::max(
                  0.0, var_mean(0, 0) + var_mean(1, 1) + var_mean(2, 2)));
    r.sigma_D = sigma_tr_c / 3.0;
    const double D_excess =
        -1.0 / 3.0 - r.mc_D;  // positive when D < -1/3 (entangled)
    r.significance_D =
        (r.sigma_D > 0.0 && D_excess > 0.0) ? D_excess / r.sigma_D : 0.0;

    r.mc_m12 = m12FromCij(r.mc_cij);
    // sigma^2[m12] = sum[i,j] (d m12 / d C_ij)^2 sigma^2[C_ij]
    //
    // Conservative rough bound: replaces each (d m12 / d C_ij)^2 by 2,
    // valid as an order-of-magnitude estimate near the Bell boundary
    // (m12 ~ 1). Exact propagation requires eigenvectors of C*C^T.
    r.sigma_m12 = std::sqrt(2.0 * r.sigma_cij.cwiseProduct(r.sigma_cij).sum());
    r.significance_bell = (r.sigma_m12 > 0.0 && r.mc_m12 > 1.0)
                              ? (r.mc_m12 - 1.0) / r.sigma_m12
                              : 0.0;

    return r;
}

// -----------------------------------------------------------------------
// computeLumiScan:
//
// The MC was run with N_MC = cfg.n_events (unweighted accepted events).
// The physical event count at luminosity L [ab^-1] is:
//
//   N(L) = sigma_eff [fb] * L [fb^-1]
//        = sigma_eff [fb] * L [ab^-1] * 1e3   (1 ab^-1 = 1000 fb^-1)
//
// Statistical uncertainties scale as 1/sqrt(N), so:
//
//   sigma_X(L)      = sigma_X(N_MC) * sqrt(N_MC / N(L))
//   significance(L) = significance(N_MC) * sqrt(N(L) / N_MC)
//                   = significance(N_MC) * sqrt(sigma_eff * L * 1e3 / N_MC)
// -----------------------------------------------------------------------
std::vector<LumiScanPoint> computeLumiScan(const MCConfig &cfg,
                                           double sigma_eff_fb,
                                           long long n_accepted,
                                           const ReconstructedMC &r) {
    std::vector<LumiScanPoint> lumi_scan;
    if (!(cfg.L_scan_min_ab < cfg.L_scan_max_ab) || cfg.n_L_points <= 1) {
        return lumi_scan;
    }

    const double N_MC = static_cast<double>(n_accepted);
    const double dL = (cfg.L_scan_max_ab - cfg.L_scan_min_ab) /
                      static_cast<double>(cfg.n_L_points - 1);

    std::cout << std::format(
        "\n-- luminosity scan [{:.3f}, {:.3f}] ab^-1, {} points\n",
        cfg.L_scan_min_ab, cfg.L_scan_max_ab, cfg.n_L_points);

    lumi_scan.reserve(cfg.n_L_points);
    for (int p = 0; p < cfg.n_L_points; ++p) {
        const double L_ab = cfg.L_scan_min_ab + p * dL;
        const double L_fb = L_ab * 1.0e3;        // ab^-1 -> fb^-1
        const double N_L = sigma_eff_fb * L_fb;  // physical event count

        if (N_L <= 0.0) {
            lumi_scan.push_back({L_ab, 0.0, 0.0, 0.0});
            continue;
        }

        const double scale = std::sqrt(N_L / N_MC);
        const double sig_C_L = r.significance_concurrence * scale;
        const double sig_D_L = r.significance_D * scale;
        const double sig_bell_L = r.significance_bell * scale;
        lumi_scan.push_back({L_ab, sig_C_L, sig_D_L, sig_bell_L});
    }

    return lumi_scan;
}
}  // namespace gagatt
