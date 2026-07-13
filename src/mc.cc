#include "mc.h"
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <random>
#include <vector>
#include "constants.h"
#include "photon.h"
#include "spin_density.h"

namespace gagatt {
// ----------------------------------------------------------------------
//  Kinematic prefactor:  (beta Nc / 32π s_hat) |A_C|²
//  |A_C|^2 = (COUPLING_FACTOR)² / (1 - beta^2 cos_th^2)^2
//  is already built into norm_factor through polCoeffsForHelicity (see
//  amplitude.cc)
//
//  partialXsec returns  d sigma_hat_{l1 l2} / d cos_th  summed over all
//  helicities with weights w^{l1 l2}:
//
//    (beta Nc / 32 pi s_hat) |A_C|² (C_1^w + C_3^w)
// ----------------------------------------------------------------------
double partialXsec(double sqrt_s_hat, double cos_th,
                   const SDMatrixCoefficients &sdc) {
    if (sdc.norm_factor <= 0.0) { return 0.0; }

    const double s_hat = sqrt_s_hat * sqrt_s_hat;
    const double r = 4.0 * MTOP2 / s_hat;
    if (r >= 1.0) { return 0.0; }

    const double beta = std::sqrt(1.0 - r);
    const double denom = 1.0 - beta * beta * cos_th * cos_th;
    if (denom <= 0.0) { return 0.0; }

    // prefactor: beta Nc / (32 pi s_hat)
    const double prefac = beta * NC / (32.0 * std::numbers::pi * s_hat);

    // sdc.norm_factor = C_1^w + C_3^w (already includes |A_C|^2 via
    // overall_fac^2) BUT overall_fac is COUPLING_FACTOR/denom, so norm_factor
    // already has the factor |A_C|² baked in. We must NOT multiply by AC2
    // again.
    return prefac * sdc.norm_factor;
}

// -----------------------------------------------------------------------
// eventRate: d^2 N / (d sqrt_s_hat  d cos_th)
//
//   = (d sigma_hat / d cos_th) [GeV^{-2}]
//   * L_tot(z) [GeV^{-2}]        (photon-photon luminosity, per unit z)
//   * L_ee [fb^{-1}]
//   * GEV2_TO_FB [fb / GeV^{-2}]
//
// The factor L_tot already carries the 1/sqrt_s = 1/(sqrt_s * z) Jacobian
// from d tau --> d sqrt_tau = d z, applied inside computePolCor (* 2z).
// An additional 1/sqrt_s factor converts dL/dz (per unit z=sqrt_tau) to
// dL/d(sqrt_s_hat): since z = sqrt_s_hat / sqrt_s, dz = d(sqrt_s_hat)/sqrt_s,
// so dN/d(sqrt_s_hat) = dN/dz * (1/sqrt_s).
// -----------------------------------------------------------------------
static double eventRate(double sqrt_s_hat, double cos_th,
                        const SDMatrixCoefficients &sdc,
                        double L_tot,      // [GeV^{-2}], per unit z
                        double sqrt_s,     // [GeV]
                        double L_ee_fb) {  // [fb^{-1}]
    const double xsec = partialXsec(sqrt_s_hat, cos_th, sdc);  // [GeV^{-2}]
    if (xsec <= 0.0 || L_tot <= 0.0) { return 0.0; }

    // dN / dz d cos_th = xsec [GeV^{-2}] * L_tot [GeV^{-2}] * L_ee [fb^{-1}]
    //                    * GEV2_TO_FB [fb/GeV^{-2}]
    // dN / d(sqrt_s_hat) d cos_th = (dN/dz d cos_th) / sqrt_s
    return xsec * L_tot * L_ee_fb * GEV2_TO_FB / sqrt_s;
}

// ----------------------------------------------------------------------
// sampleDecayAngles:
//   Sample (n_+, n_-) -- unit vectors for the two decay leptons in
//   their respective top / anti-top helicity frames -- from the joint
//   angular distribution:
//
//     P(n+, n-) = (1/4pi)^2 * [1 + B+.n+ + B-.n- + C_ij n+^i n-^j]
//
//   using accept/reject with envelope w_max.
//   Returns cos(phi) = n+ . n-.
// ----------------------------------------------------------------------
double sampleDecayAngles(const SDMatrixCoefficients &sdc,
                         std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    std::uniform_real_distribution<double> phi_dist(0.0,
                                                    2.0 * std::numbers::pi);
    std::uniform_real_distribution<double> uni01(0.0, 1.0);

    // Conservative envelope: w_max = 1 + |B+| + |B-| + ||C||_F
    const double w_max = 1.0 + sdc.bp.norm() + sdc.bm.norm() + sdc.cc.norm();

    for (;;) {
        // Sample n+ uniformly on S^2
        const double cos_th_p = uni(rng);
        const double sin_th_p = std::sqrt(1.0 - cos_th_p * cos_th_p);
        const double phi_p = phi_dist(rng);
        const Eigen::Vector3d np(sin_th_p * std::cos(phi_p),
                                 sin_th_p * std::sin(phi_p), cos_th_p);

        // Sample n- uniformly on S^2
        const double cos_th_m = uni(rng);
        const double sin_th_m = std::sqrt(1.0 - cos_th_m * cos_th_m);
        const double phi_m = phi_dist(rng);
        const Eigen::Vector3d nm(sin_th_m * std::cos(phi_m),
                                 sin_th_m * std::sin(phi_m), cos_th_m);

        // Joint weight (always >= 0 when w_max is correct)
        const double w =
            1.0 + sdc.bp.dot(np) + sdc.bm.dot(nm) + np.dot(sdc.cc * nm);

        // Accept/reject
        if (uni01(rng) * w_max <= w) {
            return np.dot(nm);  // cos(phi)
        }
    }
}

// ----------------------------------------------------------------------
// Main MC runner
// ----------------------------------------------------------------------
MCResult runMC(const MCConfig &cfg) {
    const double sqrts_max = (cfg.sqrts_max > 0.0) ? cfg.sqrts_max : cfg.sqrt_s;
    const double sqrts_min = cfg.sqrts_min;

    const double pc1 = -cfg.pe1;  // PePc = -1
    const double pc2 = -cfg.pe2;

    std::cout << std::format(
        "gagatt_mc: sqrt_s={:.0f} GeV, pe1={:+.2f}, pe2={:+.2f}, "
        "x={:.1f}, L={:.0f} fb^-1\n",
        cfg.sqrt_s, cfg.pe1, cfg.pe2, cfg.x, cfg.L_ee_fb);
    std::cout << std::format(
        "           sqrt_s_hat in [{:.1f}, {:.1f}] GeV, "
        "cos_th in [{:.2f}, {:.2f}]\n",
        sqrts_min, sqrts_max, cfg.cos_th_min, cfg.cos_th_max);

    // ------------------------------------------------------------------
    // Phase 1: precompute luminosity-weight cache
    // ------------------------------------------------------------------
    const double d_sqrts =
        (sqrts_max - sqrts_min) / static_cast<double>(cfg.n_sqrts);
    const double d_cos =
        (cfg.cos_th_max - cfg.cos_th_min) / static_cast<double>(cfg.n_cos);

    struct CacheEntry {
        LumiWeights lw;
        double L_tot;
    };
    std::vector<CacheEntry> cache(cfg.n_sqrts);

    std::cout << "-- precomputing lumi cache ...\n";
    for (int j = 0; j < cfg.n_sqrts; ++j) {
        const double ssh = sqrts_min + (j + 0.5) * d_sqrts;
        const double z = ssh / cfg.sqrt_s;
        const auto [lw, L_tot] =
            lumiWeightsAndTotal(z, cfg.x, cfg.pe1, pc1, cfg.pe2, pc2);
        cache[j].lw = lw;
        cache[j].L_tot = L_tot;

        if ((j + 1) % 500 == 0) {
            std::cout << std::format("  lumi cache: {}/{}\n", j + 1,
                                     cfg.n_sqrts);
        }
    }

    // ------------------------------------------------------------------
    // Phase 2: build 2-D weight table (used for CDF importance sampling)
    //   weight[i*n_sqrts+j] = dN / (d sqrt_s_hat d cos_th) * d_sqrts * d_cos
    //                       = number of expected events in bin (i,j)
    // ------------------------------------------------------------------
    const int N = cfg.n_cos * cfg.n_sqrts;
    std::vector<double> weights(N, 0.0);

    // Per-bin theory quantities (for weighted averages)
    std::vector<double> theory_neg(N, 0.0);
    std::vector<double> theory_con(N, 0.0);
    std::vector<double> theory_trc(N, 0.0);  // Tr[C] for <cos phi>

    // Store sdc per bin for decay sampling (reused in Phase 4)
    std::vector<SDMatrixCoefficients> sdc_cache;
    sdc_cache.reserve(N);

    std::cout << "-- building weight table ...\n";
    double total_weight = 0.0;
    for (int i = 0; i < cfg.n_cos; ++i) {
        const double cos_th = cfg.cos_th_min + (i + 0.5) * d_cos;
        for (int j = 0; j < cfg.n_sqrts; ++j) {
            const double ssh = sqrts_min + (j + 0.5) * d_sqrts;
            const auto &ce = cache[j];
            const auto sdc = SDMatrixCoefficients(ssh, cos_th, ce.lw);
            const auto rho = spinDensityMatrix(sdc);

            const double rate =
                eventRate(ssh, cos_th, sdc, ce.L_tot, cfg.sqrt_s, cfg.L_ee_fb) *
                d_sqrts * d_cos;

            const int idx = i * cfg.n_sqrts + j;
            weights[idx] = std::max(0.0, rate);
            theory_neg[idx] = negativity(rho);
            theory_con[idx] = getConcurrence(rho);
            theory_trc[idx] = sdc.cc.trace();  // Tr[C]
            sdc_cache.push_back(sdc);

            total_weight += weights[idx];
        }
        if ((i + 1) % 20 == 0) {
            std::cout << std::format("  weight table: {}/{}\n", i + 1,
                                     cfg.n_cos);
        }
    }

    if (total_weight <= 0.0) {
        std::cerr
            << "ERROR: total weight is zero — no events above threshold.\n";
        return {};
    }

    // Luminosity-weighted theory predictions
    double tw_neg = 0.0, tw_con = 0.0, tw_trc = 0.0;
    for (int k = 0; k < N; ++k) {
        tw_neg += weights[k] * theory_neg[k];
        tw_con += weights[k] * theory_con[k];
        tw_trc += weights[k] * theory_trc[k];
    }
    const double theory_cos_phi = -(tw_trc / total_weight) / 9.0;

    std::cout << std::format("-- total expected events: {:.3e}\n",
                             total_weight);
    std::cout << std::format("-- theory <cos phi>: {:.6f}\n", theory_cos_phi);
    std::cout << std::format("-- theory <negativity>:   {:.6f}\n",
                             tw_neg / total_weight);
    std::cout << std::format("-- theory <concurrence>:  {:.6f}\n",
                             tw_con / total_weight);

    // ------------------------------------------------------------------
    // Phase 3: set up discrete distribution for bin sampling (CDF)
    // ------------------------------------------------------------------
    // std::discrete_distribution samples index k with probability
    // proportional to weights[k].
    std::mt19937_64 rng(cfg.seed == 0 ? std::random_device{}()
                                      : static_cast<uint64_t>(cfg.seed));

    std::discrete_distribution<int> bin_dist(weights.begin(), weights.end());

    // ------------------------------------------------------------------
    // Phase 4: event loop — sample bin, then sample decay angles
    // ------------------------------------------------------------------
    std::cout << std::format("-- running {} MC events ...\n", cfg.n_events);

    double sum_cos_phi = 0.0;
    // double sum_cos_phi2 = 0.0;
    long long n_accepted = 0;

    const long long print_every = std::max(1LL, cfg.n_events / 10);

    while (n_accepted < cfg.n_events) {
        // Draw a bin index proportional to its expected event count
        const int k = bin_dist(rng);

        // sdc for this bin: already computed in Phase 2
        const SDMatrixCoefficients &sdc = sdc_cache[k];

        // Sample decay angles from the joint distribution and get cos(phi)
        const double cos_phi = sampleDecayAngles(sdc, rng);

        sum_cos_phi += cos_phi;
        // sum_cos_phi2 += cos_phi * cos_phi;
        ++n_accepted;

        if (n_accepted % print_every == 0) {
            std::cout << std::format("  events: {}/{}\n", n_accepted,
                                     cfg.n_events);
        }
    }

    // ------------------------------------------------------------------
    // Phase 5: aggregate statistics -> MCResult
    // ------------------------------------------------------------------
    const double n = static_cast<double>(n_accepted);
    const double mean = sum_cos_phi / n;
    // // Variance of the mean: Var(<cos phi>) = (E[x^2] - E[x]^2) / n
    // const double var_x = sum_cos_phi2 / n - mean * mean;
    // const double sigma = (n > 1.0) ? std::sqrt(std::max(0.0, var_x) / n) : 0.0;
    // // Significance: |<cos phi>| / sigma  (null: <cos phi> = 0)
    // const double signif = (sigma > 0.0) ? std::abs(mean) / sigma : 0.0;

    MCResult res;
    res.n_events_generated = n_accepted;
    res.mean_cos_phi = mean;
    // res.sigma_cos_phi = sigma;
    // res.significance = signif;
    res.theory_cos_phi = theory_cos_phi;
    res.theory_negativity = tw_neg / total_weight;
    res.theory_concurrence = tw_con / total_weight;
    res.total_xsec_fb = total_weight / cfg.L_ee_fb;

    std::cout << "\n-- MC results --\n";
    std::cout << std::format("   N events generated  : {}\n",
                             res.n_events_generated);
    std::cout << std::format("   <cos phi>  (MC)     : {:.6f} +/- {:.6f}\n",
                             res.mean_cos_phi, res.sigma_cos_phi);
    std::cout << std::format("   <cos phi>  (theory) : {:.6f}\n",
                             res.theory_cos_phi);
    std::cout << std::format("   significance        : {:.2f} sigma\n",
                             res.significance);
    std::cout << std::format("   negativity          : {:.6f}\n",
                             res.theory_negativity);
    std::cout << std::format("   concurrence         : {:.6f}\n",
                             res.theory_concurrence);
    std::cout << std::format("   total xsec          : {:.4f} fb\n",
                             res.total_xsec_fb);

    return res;
}
}  // namespace gagatt
