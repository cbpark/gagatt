#include "mc.h"
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <random>
#include <utility>
#include <vector>
#include "constants.h"
#include "photon.h"
#include "spin_density.h"

namespace gagatt {
// -----------------------------------------------------------------------
// partialXsec:
//   d sigma_hat / d cos_th  (helicity-summed, luminosity-weighted)
//   = (beta Nc / 32 pi s_hat) * |A_C|^2 * (C_1^w + C_3^w)
//
// sdc.norm_factor = C_1^w + C_3^w, with |A_C|^2 already absorbed via
// overall_fac^2 = (COUPLING_FACTOR / denom)^2 inside polCoeffsForHelicity.
// Units: [GeV^{-2}]
// -----------------------------------------------------------------------
double partialXsec(double sqrt_s_hat, double cos_th,
                   const SDMatrixCoefficients &sdc) {
    // Weighted coeffs (C_1^w, C_3^w) will be constructed by
    // `SDMatrixCoefficients::SDMatrixCoefficients(double sqrt_s_hat, double
    // cos_th, const LumiWeights &lw)` in `src/spin_density.h`
    if (sdc.norm_factor <= 0.0) { return 0.0; }

    const double s_hat = sqrt_s_hat * sqrt_s_hat;
    const double r = 4.0 * MTOP2 / s_hat;
    if (r >= 1.0) { return 0.0; }

    const double beta = std::sqrt(1.0 - r);
    const double denom = 1.0 - beta * beta * cos_th * cos_th;
    if (denom <= 0.0) { return 0.0; }

    // (beta Nc) / (32 pi s_hat)
    const double prefac = beta * NC / (32.0 * std::numbers::pi * s_hat);
    return prefac * sdc.norm_factor;  // [GeV^{-2}]
}

// -----------------------------------------------------------------------
// eventRate:
//   d^2 N / (d sqrt_s_hat  d cos_th)
//
//   = partialXsec [GeV^{-2}]
//     * L_tot(z)  [GeV^{-2}]   (photon lumi per unit z = sqrt_tau)
//     * GEV2_TO_FB
//     / sqrt_s    [GeV]         (Jacobian: dz = d(sqrt_s_hat)/sqrt_s)
// -----------------------------------------------------------------------
double eventRate(double sqrt_s_hat, double cos_th,
                 const SDMatrixCoefficients &sdc, double L_tot, double sqrt_s) {
    const double xsec = partialXsec(sqrt_s_hat, cos_th, sdc);
    if (xsec <= 0.0 || L_tot <= 0.0) { return 0.0; }

    return xsec * L_tot * GEV2_TO_FB / sqrt_s;
}

// -----------------------------------------------------------------------
// sampleDecayAngles:
//   Draw (q+, q-) — unit vectors for l+ and l- in the top / anti-top
//   rest frames — from the joint angular distribution:
//
//     P(q+, q-) = (1/4pi)^2 * [1 + B+.q+  -  B-.q-  -  q+.C.q-]
//
//   Signs:  +B+.q+  (top polarization)
//           -B-.q-  (anti-top polarization, minus sign)
//           -q+.C.q- (spin correlation, minus sign)
//
//   Accept/reject with envelope w_max = 1 + |B+| + |B-| + ||C||_F.
//   Returns the pair {q+, q-} as unit Eigen::Vector3d.
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

        // Weight:  1 + B+.q+  -  B-.q-  -  q+.C.q-
        const double w =
            1.0 + sdc.bp.dot(qp) - sdc.bm.dot(qm) - qp.dot(sdc.cc * qm);

        if (uni01(rng) * w_max <= w) { return {qp, qm}; }
    }
}

// -----------------------------------------------------------------------
// reconstructRho:
//   Build 4x4 density matrix from C_ij with B+ = B- = 0 (LO QED).
//   rho = (1/4)(I2xI2 + sum_{i,j} C_ij sigma_i x sigma_j)
// -----------------------------------------------------------------------
Matrix4cd reconstructRho(const Eigen::Matrix3d &cij) {
    using namespace Basis;
    Matrix4cd rho = I2I2;
    rho.noalias() += cij(0, 0) * S1S1 + cij(0, 1) * S1S2 + cij(0, 2) * S1S3 +
                     cij(1, 0) * S2S1 + cij(1, 1) * S2S2 + cij(1, 2) * S2S3 +
                     cij(2, 0) * S3S1 + cij(2, 1) * S3S2 + cij(2, 2) * S3S3;
    return 0.25 * rho;
}

// -----------------------------------------------------------------------
// m12FromCij:
//   Horodecki parameter m12 = two largest eigenvalues of C * C^T.
// -----------------------------------------------------------------------
double m12FromCij(const Eigen::Matrix3d &cij) {
    const Eigen::Matrix3d M = cij * cij.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(
        M, Eigen::EigenvaluesOnly);
    const auto ev = solver.eigenvalues();  // ascending order
    return ev(2) + ev(1);
}

// -----------------------------------------------------------------------
// Main MC runner
// -----------------------------------------------------------------------
MCResult runMC(const MCConfig &cfg) {
    const double sqrts_max = (cfg.sqrts_max > 0.0) ? cfg.sqrts_max : cfg.sqrt_s;
    const double sqrts_min = cfg.sqrts_min;  // ~ 2 M_t

    const double pc1 = -cfg.pe1;  // PePc = -1
    const double pc2 = -cfg.pe2;

    std::cout << std::format(
        "gagatt_mc: sqrt_s={:.0f} GeV, pe1={:+.2f}, pe2={:+.2f}, x={:.1f}\n",
        cfg.sqrt_s, cfg.pe1, cfg.pe2, cfg.x);
    std::cout << std::format(
        "  sqrt_s_hat in [{:.1f}, {:.1f}] GeV, "
        "cos_th in [{:.2f}, {:.2f}]\n",
        sqrts_min, sqrts_max, cfg.cos_th_min, cfg.cos_th_max);

    // ------------------------------------------------------------------
    // Phase 1: z-cache — precompute (LumiWeights, L_tot) per sqrt_s_hat bin
    // ------------------------------------------------------------------
    const double d_sqrts =
        (sqrts_max - sqrts_min) / static_cast<double>(cfg.n_sqrts);
    const double d_cos =
        (cfg.cos_th_max - cfg.cos_th_min) / static_cast<double>(cfg.n_cos);

    struct ZCacheEntry {
        LumiWeights lw;
        double L_tot;
    };
    std::vector<ZCacheEntry> zcache(cfg.n_sqrts);

    std::cout << "-- precomputing lumi cache ...\n";
    for (int j = 0; j < cfg.n_sqrts; ++j) {
        const double sqrt_s_hat = sqrts_min + (j + 0.5) * d_sqrts;
        const double z = sqrt_s_hat / cfg.sqrt_s;
        const auto [lw, L_tot] =
            lumiWeightsAndTotal(z, cfg.x, cfg.pe1, pc1, cfg.pe2, pc2);
        zcache[j] = {lw, L_tot};
        if ((j + 1) % 500 == 0)
            std::cout << std::format("  lumi cache: {}/{}\n", j + 1,
                                     cfg.n_sqrts);
    }

    // ------------------------------------------------------------------
    // Phase 2: build 2-D weight table; cache sdc; accumulate means
    // ------------------------------------------------------------------
    const int N_bins = cfg.n_cos * cfg.n_sqrts;
    std::vector<double> bin_weights(N_bins, 0.0);
    std::vector<SDMatrixCoefficients> sdc_cache;
    sdc_cache.reserve(N_bins);

    double tw_neg = 0.0, tw_con = 0.0, tw_trc = 0.0, tw_m12 = 0.0;
    double total_weight = 0.0;

    std::cout << "-- building weight table ...\n";
    for (int i = 0; i < cfg.n_cos; ++i) {
        const double cos_th = cfg.cos_th_min + (i + 0.5) * d_cos;
        for (int j = 0; j < cfg.n_sqrts; ++j) {
            const double sqrt_s_hat = sqrts_min + (j + 0.5) * d_sqrts;
            const SDMatrixCoefficients sdc(sqrt_s_hat, cos_th, zcache[j].lw);
            const Matrix4cd rho = spinDensityMatrix(sdc);

            const double rate =
                eventRate(sqrt_s_hat, cos_th, sdc, zcache[j].L_tot, cfg.sqrt_s) *
                d_sqrts * d_cos;

            const int idx = i * cfg.n_sqrts + j;
            bin_weights[idx] = std::max(0.0, rate);
            sdc_cache.push_back(sdc);

            total_weight += bin_weights[idx];
            tw_neg += bin_weights[idx] * negativity(rho);
            tw_con += bin_weights[idx] * getConcurrence(rho);
            tw_trc += bin_weights[idx] * sdc.cc.trace();
            tw_m12 += bin_weights[idx] * horodeckiMeasure(sdc);
        }
        if ((i + 1) % 20 == 0)
            std::cout << std::format("  weight table: {}/{}\n", i + 1,
                                     cfg.n_cos);
    }

    if (total_weight <= 0.0) {
        std::cerr << "ERROR: total weight is zero — check kinematics.\n";
        return {};
    }

    const double theory_tr_c = tw_trc / total_weight;
    const double theory_neg = tw_neg / total_weight;
    const double theory_con = tw_con / total_weight;
    const double theory_m12 = tw_m12 / total_weight;

    std::cout << std::format("-- total expected events : {:.3e}\n",
                             total_weight);
    std::cout << std::format("-- theory Tr[C]          : {:+.6f}\n",
                             theory_tr_c);
    std::cout << std::format("-- theory <cos phi>      : {:+.6f}\n",
                             -theory_tr_c / 9.0);
    std::cout << std::format("-- theory negativity     : {:.6f}\n", theory_neg);
    std::cout << std::format("-- theory concurrence    : {:.6f}\n", theory_con);
    std::cout << std::format("-- theory m12            : {:.6f}\n", theory_m12);

    // ------------------------------------------------------------------
    // Phase 3: discrete sampler proportional to expected event counts
    // ------------------------------------------------------------------
    std::mt19937_64 rng(cfg.seed == 0 ? std::random_device{}()
                                      : static_cast<uint64_t>(cfg.seed));
    std::discrete_distribution<int> bin_dist(bin_weights.begin(),
                                             bin_weights.end());

    // ------------------------------------------------------------------
    // Phase 4: event loop
    //
    //   For each event, draw a (cos_th, sqrt_s_hat) bin, then sample
    //   decay directions (q+, q-).
    //
    //   Accumulate:
    //     S_ij      = sum  q+_i * q-_j          (first moment)
    //     S2_ij     = sum (q+_i * q-_j)^2       (second moment, for variance)
    //
    //   After N events:
    //     <q+_i q-_j>  = S_ij / N
    //     C_ij^MC      = -9 * S_ij / N
    //     Var[<q+_i q-_j>] = (S2_ij/N - (S_ij/N)^2) / N
    //     sigma[C_ij]  = 9 * sqrt(Var[<q+_i q-_j>])
    // ------------------------------------------------------------------
    std::cout << std::format("-- running {} MC events ...\n", cfg.n_events);

    Eigen::Matrix3d S1_qpqm = Eigen::Matrix3d::Zero();  // sum q+_i q-_j
    Eigen::Matrix3d S2_qpqm = Eigen::Matrix3d::Zero();  // sum (q+_i q-_j)^2
    long long n_accepted = 0;

    const long long print_every = std::max(1LL, cfg.n_events / 10);

    while (n_accepted < cfg.n_events) {
        const int k = bin_dist(rng);
        const SDMatrixCoefficients &sdc = sdc_cache[k];

        // Sample q+ and q- from the joint distribution
        const auto [qp, qm] = sampleDecayAngles(sdc, rng);

        // Accumulate outer product and its element-wise square
        const Eigen::Matrix3d outer = qp * qm.transpose();
        S1_qpqm += outer;
        S2_qpqm += outer.cwiseProduct(outer);

        ++n_accepted;
        if (n_accepted % print_every == 0)
            std::cout << std::format("  events: {}/{}\n", n_accepted,
                                     cfg.n_events);
    }

    // ------------------------------------------------------------------
    // Phase 5: reconstruct C_ij and all derived quantities
    // ------------------------------------------------------------------
    const double n = static_cast<double>(n_accepted);

    // <q+_i q-_j>
    const Eigen::Matrix3d mean_qpqm = S1_qpqm / n;

    // C_ij^MC = -9 * <q+_i q-_j>
    const Eigen::Matrix3d mc_cij = -9.0 * mean_qpqm;

    // Variance of the mean per element: Var[<q+_i q-_j>] = (E[x^2] - E[x]^2) /
    // N
    const Eigen::Matrix3d mean_qpqm2 = S2_qpqm / n;
    const Eigen::Matrix3d var_mean =
        (mean_qpqm2 - mean_qpqm.cwiseProduct(mean_qpqm)) / n;

    // sigma[C_ij] = 9 * sqrt(Var[<q+_i q-_j>])
    const Eigen::Matrix3d sigma_cij = 9.0 * var_mean.cwiseSqrt().cwiseAbs();

    // Tr[C] and its uncertainty:
    //   sigma^2[Tr[C]] = 81 * (Var[<q+_n q-_n>] + Var[<q+_r q-_r>] + Var[<q+_k
    //   q-_k>])
    const double mc_tr_c = mc_cij.trace();
    const double sigma_tr_c =
        9.0 * std::sqrt(std::max(
                  0.0, var_mean(0, 0) + var_mean(1, 1) + var_mean(2, 2)));

    // D = Tr[C]/3 ;  significance vs null (D = -1/3)
    const double sigma_D = sigma_tr_c / 3.0;
    const double D_val = mc_tr_c / 3.0;
    const double D_excess = -D_val - 1.0 / 3.0;  // positive when D < -1/3
    const double significance_D =
        (sigma_D > 0.0 && D_excess > 0.0) ? D_excess / sigma_D : 0.0;
    std::cout << std::format(" D_val              : {:+.6f}\n", D_val);
    std::cout << std::format(" D_excess (|D|-1/3) : {:+.6f}\n", D_excess);
    std::cout << std::format(" sigma_D (N_MC)     : {:.6f}\n", sigma_D);
    std::cout << std::format(
        " sig_D at N_MC      : {:.2f} sigma\n",
        (sigma_D > 0.0 && D_excess > 0.0) ? D_excess / sigma_D : 0.0);

    // Reconstruct density matrix from C_ij^MC (B+ = B- = 0 at LO)
    const Matrix4cd mc_rho = reconstructRho(mc_cij);
    const double mc_negativity = negativity(mc_rho);
    const double mc_concurrence = getConcurrence(mc_rho);
    const double mc_m12 = m12FromCij(mc_cij);

    // sigma^2[m12] = sum[i,j] (d m12 / d C_ij)^2 sigma^2[C_ij]
    //
    // Conservative rough bound: replaces each (d m12 / d C_ij)^2 by 2,
    // valid as an order-of-magnitude estimate near the Bell boundary (m12 ~ 1).
    // Exact propagation requires eigenvectors of C*C^T.
    const double sigma_m12 =
        std::sqrt(2.0 * sigma_cij.cwiseProduct(sigma_cij).sum());
    const double significance_bell =
        (sigma_m12 > 0.0 && mc_m12 > 1.0) ? (mc_m12 - 1.0) / sigma_m12 : 0.0;

    // ------------------------------------------------------------------
    // Phase 6: print results
    // ------------------------------------------------------------------
    std::cout << "\n-- MC results --\n";
    std::cout << std::format("  N events generated   : {}\n", n_accepted);

    const double sigma_prod_fb = total_weight;
    const double sigma_eff_fb = sigma_prod_fb * BRLL;
    std::cout << std::format(" production xsec (ee->gaga->tt): {:.4f} fb\n",
                             sigma_prod_fb);
    std::cout << std::format(" effective xsec  (* BR_ll)     : {:.4f} fb\n",
                             sigma_eff_fb);
    std::cout << std::format(" BR(tt->ll)                    : {:.4f}\n", BRLL);
    std::cout << std::format(" N(L=0.01 ab^-1)    : {:.1f} events\n",
                             sigma_eff_fb * 10.0);

    std::cout << "\n  Reconstructed C_ij  (rows: n,r,k; cols: n,r,k)\n";
    const std::array<const char *, 3> ax = {"n", "r", "k"};
    for (int a = 0; a < 3; ++a) {
        std::cout << std::format("    C_{0}j :", ax[a]);
        for (int b = 0; b < 3; ++b)
            std::cout << std::format("  {:+.4f}+/-{:.4f}", mc_cij(a, b),
                                     sigma_cij(a, b));
        std::cout << "\n";
    }

    std::cout << std::format("\n  Tr[C]  (MC)          : {:+.6f} +/- {:.6f}\n",
                             mc_tr_c, sigma_tr_c);
    std::cout << std::format("  Tr[C]  (theory)      : {:+.6f}\n", theory_tr_c);
    std::cout << std::format("  D=Tr[C]/3  (MC)      : {:+.6f} +/- {:.6f}\n",
                             mc_tr_c / 3.0, sigma_D);
    std::cout << std::format("  significance(D)      : {:.2f} sigma\n",
                             significance_D);

    std::cout << std::format("\n  negativity  (MC)     : {:.6f}\n",
                             mc_negativity);
    std::cout << std::format("  negativity  (theory) : {:.6f}\n", theory_neg);
    std::cout << std::format("  concurrence (MC)     : {:.6f}\n",
                             mc_concurrence);
    std::cout << std::format("  concurrence (theory) : {:.6f}\n", theory_con);

    // std::cout << std::format("\n  m12         (MC)     : {:.6f}\n", mc_m12);
    // std::cout << std::format("  m12         (theory) : {:.6f}\n",
    // theory_m12); std::cout << std::format("  significance(Bell)   : {:.2f}
    // sigma\n",
    //                          significance_bell);
    std::cout << std::format(" m12 (MC)          : {:.6f}\n", mc_m12);
    std::cout << std::format(" m12 - 1           : {:.6f}\n", mc_m12 - 1.0);
    std::cout << std::format(" sigma_m12 (N_MC)  : {:.6f}\n", sigma_m12);
    std::cout << std::format(
        " sig_Bell at N_MC  : {:.2f} sigma\n",
        (mc_m12 > 1.0 && sigma_m12 > 0.0) ? (mc_m12 - 1.0) / sigma_m12 : 0.0);

    // ------------------------------------------------------------------
    // Phase 7: luminosity scan
    //
    // The MC was run with N_MC = cfg.n_events (unweighted accepted events).
    // The physical event count at luminosity L [ab^-1] is:
    //
    //   N(L) = sigma_eff [fb] * L [fb^-1]
    //         = sigma_eff [fb] * L [ab^-1] * 1e3  (1 ab^-1 = 1000 fb^-1)
    //
    // Statistical uncertainties scale as 1/sqrt(N), so:
    //
    //   sigma_X(L)     = sigma_X(N_MC) * sqrt(N_MC / N(L))
    //   significance(L) = significance(N_MC) * sqrt(N(L) / N_MC)
    //                   = significance(N_MC) * sqrt(sigma_eff * L * 1e3 / N_MC)
    // ------------------------------------------------------------------
    std::vector<LumiScanPoint> lumi_scan;

    if (cfg.L_scan_min_ab < cfg.L_scan_max_ab && cfg.n_L_points > 1) {
        const double N_MC = static_cast<double>(n_accepted);

        const double dL = (cfg.L_scan_max_ab - cfg.L_scan_min_ab) /
                          static_cast<double>(cfg.n_L_points - 1);

        std::cout << std::format(
            "\n-- luminosity scan [{:.3f}, {:.3f}] ab^-1, {} points\n",
            cfg.L_scan_min_ab, cfg.L_scan_max_ab, cfg.n_L_points);
        std::cout << std::format("   N_MC = {:.0f}\n", N_MC);

        lumi_scan.reserve(cfg.n_L_points);

        for (int p = 0; p < cfg.n_L_points; ++p) {
            const double L_ab = cfg.L_scan_min_ab + p * dL;
            const double L_fb = L_ab * 1.0e3;        // ab^-1 -> fb^-1
            const double N_L = sigma_eff_fb * L_fb;  // physical event count

            if (N_L <= 0.0) {
                lumi_scan.push_back({L_ab, 0.0, 0.0});
                continue;
            }

            const double scale = std::sqrt(N_L / N_MC);

            // const double sig_D_L = significance_D * scale;
            const double sig_D_L =
                (significance_D > 0.0) ? significance_D * scale : 0.0;

            const double sig_bell_L = (mc_m12 > 1.0 && sigma_m12 > 0.0)
                                          ? (mc_m12 - 1.0) / (sigma_m12 / scale)
                                          : 0.0;

            lumi_scan.push_back({L_ab, sig_D_L, sig_bell_L});
        }

        // Print table to stdout
        std::cout << std::format(" {:>12s}  {:>16s}  {:>16s}\n", "L[ab^-1]",
                                 "sig_D[sigma]", "sig_Bell[sigma]");
        for (const auto &pt : lumi_scan) {
            std::cout << std::format(" {:12.6f}  {:16.6f}  {:16.6f}\n", pt.L_ab,
                                     pt.significance_D, pt.significance_bell);
        }
    }

    // ------------------------------------------------------------------
    // Phase 8: fill and return MCResult
    // ------------------------------------------------------------------
    MCResult res;
    res.n_events_generated = n_accepted;
    res.mc_cij = mc_cij;
    res.sigma_cij = sigma_cij;
    res.mc_tr_c = mc_tr_c;
    res.sigma_tr_c = sigma_tr_c;
    res.significance_D = significance_D;
    res.mc_concurrence = mc_concurrence;
    res.mc_negativity = mc_negativity;
    res.mc_m12 = mc_m12;
    res.significance_bell = significance_bell;
    res.theory_tr_c = theory_tr_c;
    res.theory_concurrence = theory_con;
    res.theory_negativity = theory_neg;
    res.theory_m12 = theory_m12;
    res.total_xsec_fb = total_weight;
    res.lumi_scan = std::move(lumi_scan);
    return res;
}

}  // namespace gagatt
