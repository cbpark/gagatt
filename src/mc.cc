#include "mc.h"
#include <format>
#include <iostream>
#include <random>
#include "mc_helper.h"
#include "spin_density.h"

namespace gagatt {
// -----------------------------------------------------------------------
// Main MC runner
//
//   Phase 1: buildLumiCache: z-cache of (LumiWeights, L_tot)
//   Phase 2: buildWeightTable: 2-D weight table + theory averages
//   Phase 3: runEventLoop: accept/reject sampling of (q+, q-)
//   Phase 4: print results
//   Phase 5: computeLumiScan: significance vs. integrated luminosity
//   Phase 6: fill and return MCResult
// -----------------------------------------------------------------------
MCResult runMC(const MCConfig &cfg) {
    const double sqrts_max = (cfg.sqrts_max > 0.0) ? cfg.sqrts_max : cfg.sqrt_s;
    const double sqrts_min = cfg.sqrts_min;  // ~ 2 M_t

    const double pc1 = -cfg.pe1;  // pc = -pe
    const double pc2 = -cfg.pe2;

    std::cout << std::format(
        "gagatt_mc: sqrt_s={:.0f} GeV, pe1={:+.2f}, pe2={:+.2f}, x={:.1f}\n",
        cfg.sqrt_s, cfg.pe1, cfg.pe2, cfg.x);
    std::cout << std::format(
        "  sqrt_s_hat in [{:.1f}, {:.1f}] GeV, "
        "cos_th in [{:.2f}, {:.2f}]\n",
        sqrts_min, sqrts_max, cfg.cos_th_min, cfg.cos_th_max);

    const double d_sqrts =
        (sqrts_max - sqrts_min) / static_cast<double>(cfg.n_sqrts);
    const double d_cos =
        (cfg.cos_th_max - cfg.cos_th_min) / static_cast<double>(cfg.n_cos);

    // Phase 1: lumi cache
    const std::vector<ZCacheEntry> zcache =
        buildLumiCache(cfg, pc1, pc2, sqrts_min, d_sqrts);

    // Phase 2: weight table + theory averages
    const WeightTable wt =
        buildWeightTable(cfg, zcache, sqrts_min, d_sqrts, d_cos);
    if (wt.bin_weights.empty()) { return {}; }  // total_weight <= 0

    std::mt19937_64 rng(cfg.seed == 0 ? std::random_device{}()
                                      : static_cast<uint64_t>(cfg.seed));

    // Phase 3: event loop
    const EventLoopResult ev = runEventLoop(cfg.n_events, wt, rng);
    const ReconstructedMC r = reconstructFromMoments(ev);

    // ------------------------------------------------------------------
    // Phase 4: print results
    // ------------------------------------------------------------------
    std::cout << "\n-- MC results --\n";
    std::cout << std::format(" N events generated            : {}\n",
                             ev.n_accepted);

    const double sigma_prod_fb = wt.total_weight;
    const double sigma_eff_fb = sigma_prod_fb * (BRLL + BRLJ);
    std::cout << std::format(" production xsec (ee->gaga->tt)   : {:.4f} fb\n",
                             sigma_prod_fb);
    std::cout << std::format(" effective xsec (* (BR_ll + BRlj)): {:.4f} fb\n",
                             sigma_eff_fb);
    std::cout << std::format(" BR(tt->ll/lj)                    : {:.4f}\n", BRLL + BRLJ);

    std::cout << std::format("\n Concurrence (theory)  : {:.6f}\n",
                             wt.theory_concurrence);
    std::cout << std::format(" Concurrence (MC, N={}): {:.6f} +/- {:.6f}\n",
                             cfg.n_events, r.mc_concurrence,
                             r.sigma_concurrence);
    std::cout << std::format(" D           (theory)  : {:.6f}\n", wt.theory_D);
    std::cout << std::format(" D           (MC, N={}): {:.6f} +/- {:.6f}\n",
                             cfg.n_events, r.mc_D, r.sigma_D);
    std::cout << std::format(" m12         (theory)  : {:.6f}\n",
                             wt.theory_m12);
    std::cout << std::format(" m12         (MC, N={}): {:.6f} +/- {:.6f}\n",
                             cfg.n_events, r.mc_m12, r.sigma_m12);

    // Phase 5: multi-seed luminosity scan
    // rng is passed through: each seed is deterministically derived
    // from the main seed but independent of the Phase 3 event loop.
    const std::vector<LumiScanPoint> lumi_scan =
        computeLumiScan(cfg, sigma_eff_fb, wt, rng);

    // ------------------------------------------------------------------
    // Phase 6: fill and return MCResult
    // ------------------------------------------------------------------
    MCResult res;
    res.n_events_generated = ev.n_accepted;

    res.mc_concurrence = r.mc_concurrence;
    res.mc_D = r.mc_D;
    res.mc_m12 = r.mc_m12;
    res.sigma_concurrence = r.sigma_concurrence;
    res.sigma_D = r.sigma_D;
    res.sigma_m12 = r.sigma_m12;

    res.theory_concurrence = wt.theory_concurrence;
    res.theory_D = wt.theory_D;
    res.theory_m12 = wt.theory_m12;

    res.total_xsec_fb = wt.total_weight;
    res.lumi_scan = std::move(lumi_scan);
    return res;
}
}  // namespace gagatt
