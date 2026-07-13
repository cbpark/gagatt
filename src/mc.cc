#include "mc.h"
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
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

double lumiTotal(double z, double x, double pe1, double pc1, double pe2,
                 double pc2) {
    const double sc1 = sigmaC(x, pe1, pc1);
    const double sc2 = sigmaC(x, pe2, pc2);
    if (sc1 <= 0.0 || sc2 <= 0.0) { return 0.0; }

    // Re-use photonLuminosity with nullopt helicities --> returns L^unp
    // which = <00>_tau / (sigma_c1 sigma_c2).
    // L^tot = 4 L^avg = 4 L^unp
    const double L_unp = photonLuminosity(z, x, pe1, pc1, pe2, pc2, {}, {});
    return 4.0 * L_unp;
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

    // Phase 1: precompute luminosity-weight cache
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
        cache[j].lw = lumiWeights(z, cfg.x, cfg.pe1, pc1, cfg.pe2, pc2);
        cache[j].L_tot = lumiTotal(z, cfg.x, cfg.pe1, pc1, cfg.pe2, pc2);

        if ((j + 1) % 500 == 0) {
            std::cout << std::format("   lumi cache: {}/{}\n", j + 1,
                                     cfg.n_sqrts);
        }
    }

    // Phase 2: build 2D weight table for importance sampling
    // w[i][j] = event rate at (cos_th_i, sqrt_s_hat_j)
    // We flatten to 1-D for CDF sampling.
    const int N = cfg.n_cos * cfg.n_sqrts;
    std::vector<double> weights(N);
    std::vector<double> theory_neg(N, 0.0);
    std::vector<double> theory_con(N, 0.0);
    std::vector<double> theory_trc(N, 0.0);  // Tr[C] for <cos phi>

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
                eventRate(ssh, cos_th, sdc, ce.L_tot, cfg.L_ee_fb) * d_sqrts *
                d_cos;  // bin weight

            const int idx = i * cfg.n_sqrts + j;
            weights[idx] = std::max(0.0, rate);
            theory_neg[idx] = negativity(rho);
            theory_con[idx] = getConcurrence(rho);
            theory_trc[idx] = sdc.cc.trace();  // Tr[C]
            total_weight += weights[idx];
        }
        if ((i + 1) % 20 == 0) {
            std::cout << std::format("   weight table: {}/{}\n", i + 1,
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
    const double theory_cos_phi = -(tw_trc / total_weight) / 9.0;  // Eq. (5.35)

    // Phase 5: compute results
    MCResult res;

    return res;
}
}  // namespace gagatt
