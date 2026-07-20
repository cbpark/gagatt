#include "mc_helper.h"
#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <numeric>
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
static double partialXsec(double sqrt_s_hat, const SDMatrixCoefficients &sdc) {
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

            const double rate =
                eventRate(sqrt_s_hat, sdc, zcache[j].L_tot, cfg.sqrt_s) *
                d_sqrts * d_cos;

            const int idx = i * cfg.n_sqrts + j;
            wt.bin_weights[idx] = std::max(0.0, rate);
            wt.sdc_cache.push_back(sdc);

            wt.total_weight += wt.bin_weights[idx];

            // theory_concurrence = <C(rho)>,
            // the phase-space average of per-bin concurrence.
            tw_con += wt.bin_weights[idx] * getConcurrence(sdc);
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
        return {};
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

    // Frobenius norm >= operator norm,
    // so this is a valid (conservative) envelope
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
EventLoopResult runEventLoop(long long n_events, const WeightTable &wt,
                             std::mt19937_64 &rng, bool verbose) {
    std::discrete_distribution<int> bin_dist(wt.bin_weights.begin(),
                                             wt.bin_weights.end());

    EventLoopResult ev;
    ev.per_bin.resize(wt.bin_weights.size());
    const long long print_every = std::max(1LL, n_events / 10);

    while (ev.n_accepted < n_events) {
        const int k = bin_dist(rng);
        const SDMatrixCoefficients &sdc = wt.sdc_cache[k];

        // Sample q+ and q- from the joint distribution
        const auto [qp, qm] = sampleDecayAngles(sdc, rng);

        // Accumulate outer product and its element-wise square (for C_ij)
        const Eigen::Matrix3d outer = qp * qm.transpose();

        // per-bin accumulation — k identifies which (cos_th, sqrt_s_hat) bin
        ev.per_bin[k].S1_qpqm += outer;
        ev.per_bin[k].S2_qpqm += outer.cwiseProduct(outer);
        ev.per_bin[k].S1_qp += qp;
        ev.per_bin[k].S1_qm += qm;
        ev.per_bin[k].n += 1;

        ++ev.n_accepted;
        if (verbose && ev.n_accepted % print_every == 0) {
            std::cout << std::format(" events: {}/{}\n", ev.n_accepted,
                                     n_events);
        }
    }
    return ev;
}

// -----------------------------------------------------------------------
// reconstructFromMoments:
// -----------------------------------------------------------------------
ReconstructedMC reconstructFromMoments(const EventLoopResult &ev) {
    ReconstructedMC r;

    // Compute <C(rho_k)>: concurrence / D / m12 per bin, then average.
    double sum_w = 0.0, sum_wC = 0.0, sum_wD = 0.0, sum_wm12 = 0.0;
    double sum_w2_varC = 0.0;    // sum n_k^2 * sigma_C_k^2
    double sum_w2_varD = 0.0;    // sum n_k^2 * sigma_D_k^2
    double sum_w2_varm12 = 0.0;  // sum n_k^2 * sigma_m12_k^2
    for (const auto &bin : ev.per_bin) {
        if (bin.n < 1) { continue; }
        const double nk = static_cast<double>(bin.n);

        const Eigen::Matrix3d cij_k = -9.0 * (bin.S1_qpqm / nk);
        const Eigen::Vector3d bp_k = 3.0 * (bin.S1_qp / nk);
        const Eigen::Vector3d bm_k = -3.0 * (bin.S1_qm / nk);

        if (cij_k.cwiseAbs().maxCoeff() > 1.0 ||
            bp_k.cwiseAbs().maxCoeff() > 1.0 ||
            bm_k.cwiseAbs().maxCoeff() > 1.0) {
            continue;
        }

        const Matrix4cd rho_k = reconstructRho(bp_k, bm_k, cij_k);

        sum_w += nk;
        sum_wC += nk * getConcurrence(rho_k);
        sum_wD += nk * entanglementMarker(cij_k);
        sum_wm12 += nk * m12FromCij(cij_k);

        const Eigen::Matrix3d mean_qpqm_k = bin.S1_qpqm / nk;
        const Eigen::Matrix3d mean_qpqm2_k = bin.S2_qpqm / nk;
        const Eigen::Matrix3d var_mean_k =
            (mean_qpqm2_k - mean_qpqm_k.cwiseProduct(mean_qpqm_k)) / nk;

        // --- sigma_concurrence per bin ---
        // sigma_cij_k(i,j) = 9 * sqrt(Var[<qp_i qm_j>_k])
        // Conservative: sigma_con_k = 0.5 * ||sigma_cij_k||_F
        const Eigen::Matrix3d sigma_cij_k =
            9.0 * var_mean_k.cwiseMax(0.0).cwiseSqrt();
        const double sigma_con_k = 0.5 * sigma_cij_k.norm();
        sum_w2_varC += nk * nk * sigma_con_k * sigma_con_k;

        // --- sigma_D per bin ---
        // D_k = (C_nn_k - |C_rr_k + C_kk_k|) / 3
        // |dD/dC_nn| = |dD/dC_rr| = |dD/dC_kk| = 1/3
        // => sigma_D_k = (1/3) * sqrt(var_nn_k + var_rr_k + var_kk_k)
        //   where var_XY_k = 81 * var_mean_k(i,i)  [factor 9^2 from C_ij =
        //   -9*<qp qm>]
        // Combining: sigma_D_k = 3 * sqrt(var_mean_k(0,0) + var_mean_k(1,1) +
        // var_mean_k(2,2))
        const double var_diag_k = std::max(
            0.0, var_mean_k(0, 0) + var_mean_k(1, 1) + var_mean_k(2, 2));
        const double sigma_D_k = 3.0 * std::sqrt(var_diag_k);
        sum_w2_varD += nk * nk * sigma_D_k * sigma_D_k;

        // --- sigma_m12 per bin ---
        // Conservative bound: sigma_m12_k = sqrt(2) * ||sigma_cij_k||_F
        //                                 = sqrt(2) * 9 * sqrt(sum var_mean_k
        //                                 entries)
        const double sigma_m12_k = std::sqrt(2.0 * sigma_cij_k.squaredNorm());
        sum_w2_varm12 += nk * nk * sigma_m12_k * sigma_m12_k;
    }
    r.mc_concurrence = (sum_w > 0.0) ? sum_wC / sum_w : 0.0;
    r.mc_D = (sum_w > 0.0) ? sum_wD / sum_w : 0.0;
    r.mc_m12 = (sum_w > 0.0) ? sum_wm12 / sum_w : 0.0;

    // sigma^2 = (1/W)^2 * sum_k n_k^2 * sigma_X_k^2
    r.sigma_concurrence = (sum_w > 0.0) ? std::sqrt(sum_w2_varC) / sum_w : 0.0;
    r.sigma_D = (sum_w > 0.0) ? std::sqrt(sum_w2_varD) / sum_w : 0.0;
    r.sigma_m12 = (sum_w > 0.0) ? std::sqrt(sum_w2_varm12) / sum_w : 0.0;

    // significance of C > 0 (entanglement detection threshold)
    r.significance_concurrence =
        (r.sigma_concurrence > 0.0 && r.mc_concurrence > 0.0)
            ? r.mc_concurrence / r.sigma_concurrence
            : 0.0;

    const double D_excess =
        -1.0 / 3.0 - r.mc_D;  // positive when D < -1/3 (entangled)
    r.significance_D =
        (r.sigma_D > 0.0 && D_excess > 0.0) ? D_excess / r.sigma_D : 0.0;

    r.significance_m12 = (r.sigma_m12 > 0.0 && r.mc_m12 > 1.0)
                             ? (r.mc_m12 - 1.0) / r.sigma_m12
                             : 0.0;

    return r;
}

// -----------------------------------------------------------------------
// computeLumiScan:
//
// For each luminosity point L:
//   1. N(L) = sigma_eff_fb * L_fb  (physical event count)
//   2. Run cfg.n_lumi_seeds independent event loops, each with N(L)
//   events, using successive RNG states drawn from the provided rng.
//   3. For each seed reconstruct mc_concurrence, sigma_concurrence, etc.
//      and compute the per-seed significance.
//   4. Report mean +/- std-dev of significance across seeds.
//
// This correctly models the finite-statistics degradation of the
// measured concurrence at low luminosity, which sqrt(N) rescaling
// cannot capture.
// -----------------------------------------------------------------------
std::vector<LumiScanPoint> computeLumiScan(const MCConfig &cfg,
                                           double sigma_eff_fb,
                                           const WeightTable &wt,
                                           std::mt19937_64 &rng) {
    std::vector<LumiScanPoint> lumi_scan;
    if (cfg.L_scan_min_ab >= cfg.L_scan_max_ab || cfg.n_L_points <= 1) {
        return lumi_scan;
    }
    if (sigma_eff_fb <= 0.0) { return lumi_scan; }

    const double dL = (cfg.L_scan_max_ab - cfg.L_scan_min_ab) /
                      static_cast<double>(cfg.n_L_points - 1);
    const int K = std::max(1, cfg.n_lumi_seeds);

    std::cout << std::format(
        "\n-- luminosity scan [{:.3f}, {:.3f}] ab^-1, {} points, {} "
        "seeds/point\n",
        cfg.L_scan_min_ab, cfg.L_scan_max_ab, cfg.n_L_points, K);
    std::cout << std::format(
        "   {:>6}  {:>10}  {:>12}  {:>12}  {:>10}  {:>10}  {:>10}\n", "pt",
        "L[ab^-1]", "N(L)", "C_mean", "sig_C", "sigma_C", "sig_D");

    lumi_scan.reserve(cfg.n_L_points);
    for (int p = 0; p < cfg.n_L_points; ++p) {
        const double L_ab = cfg.L_scan_min_ab + p * dL;
        const double L_fb = L_ab * 1.0e3;        // ab^-1 -> fb^-1
        const double N_L = sigma_eff_fb * L_fb;  // physical event count
        const long long N_events = static_cast<long long>(std::max(1.0, N_L));

        if (N_L <= 0.0) {
            lumi_scan.push_back({L_ab, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
            std::cout << std::format(
                "   {:>6}  {:>10.4f}  {:>12.1f}  {:>12}  {:>10}  {:>10}  "
                "{:>10}\n",
                p + 1, L_ab, 0.0, "N/A", "N/A", "N/A", "N/A");
            continue;
        }

        std::vector<double> sig_C(K), sig_D(K), sig_m12(K);
        std::vector<double> con(K), scon(K);  // per-seed C and sigma_C
        for (int k = 0; k < K; ++k) {
            const uint64_t seed_k = rng();  // draw from main rng
            std::mt19937_64 rng_k(seed_k);

            const EventLoopResult ev_k = runEventLoop(N_events, wt, rng_k,
                                                      /*verbose=*/false);
            const ReconstructedMC recon_k = reconstructFromMoments(ev_k);

            con[k] = recon_k.mc_concurrence;
            scon[k] = recon_k.sigma_concurrence;
            sig_C[k] = recon_k.significance_concurrence;
            sig_D[k] = recon_k.significance_D;
            sig_m12[k] = recon_k.significance_m12;
        }

        auto mean_of = [&](const std::vector<double> &v) {
            return std::accumulate(v.begin(), v.end(), 0.0) / K;
        };
        auto stddev_of = [&](const std::vector<double> &v, double mu) {
            if (K <= 1) return 0.0;
            double acc = 0.0;
            for (double x : v) acc += (x - mu) * (x - mu);
            return std::sqrt(acc / (K - 1));  // sample std-dev
        };

        const double mC = mean_of(sig_C);
        const double mD = mean_of(sig_D);
        const double mm12 = mean_of(sig_m12);
        const double mcon = mean_of(con);
        const double mscon = mean_of(scon);

        // progress line: one row per luminosity point
        std::cout << std::format(
            "   {:>6}  {:>10.4f}  {:>12.1f}  {:>12.6f}  {:>10.4f}  {:>10.6f}  "
            "{:>10.4f}\n",
            p + 1, L_ab, N_L, mcon, mC, mscon, mD);

        lumi_scan.push_back({L_ab, mC, mD, mm12, stddev_of(sig_C, mC),
                             stddev_of(sig_D, mD), stddev_of(sig_m12, mm12)});
    }
    std::cout << std::flush;
    return lumi_scan;
}
}  // namespace gagatt
