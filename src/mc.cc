#include "mc.h"
#include <format>
#include <iostream>
#include <vector>
#include "photon.h"

namespace gagatt {
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

    // Phase 5: compute results
    MCResult res;

    return res;
}
}  // namespace gagatt
