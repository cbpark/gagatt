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
//   Phase 4: reconstructFromMoments: C_ij^MC and derived quantities
//   Phase 5: print results
//   Phase 6: computeLumiScan: significance vs. integrated luminosity
//   Phase 7: fill and return MCResult
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
    if (wt.bin_weights.empty()) { return {}; }

    // Phase 3: event loop
    std::mt19937_64 rng(cfg.seed == 0 ? std::random_device{}()
                                      : static_cast<uint64_t>(cfg.seed));
    const EventLoopResult ev = runEventLoop(cfg.n_events, wt, rng);

    // Phase 4: reconstruct C_ij and all derived quantities
    const ReconstructedMC r = reconstructFromMoments(ev);

    // ------------------------------------------------------------------
    // Phase 5: print results
    // ------------------------------------------------------------------
    std::cout << "\n-- MC results --\n";
    std::cout << std::format(" N events generated            : {}\n",
                             ev.n_accepted);

    const double sigma_prod_fb = wt.total_weight;
    const double sigma_eff_fb = sigma_prod_fb * BRLL;
    std::cout << std::format(" production xsec (ee->gaga->tt): {:.4f} fb\n",
                             sigma_prod_fb);
    std::cout << std::format(" effective xsec (* BR_ll)      : {:.4f} fb\n",
                             sigma_eff_fb);
    std::cout << std::format(" BR(tt->ll)                    : {:.4f}\n", BRLL);

    std::cout << "\n Reconstructed C_ij (rows: n,r,k; cols: n,r,k)\n";
    const std::array<const char *, 3> ax = {"n", "r", "k"};
    for (int a = 0; a < 3; ++a) {
        std::cout << std::format(" C_{0}j :", ax[a]);
        for (int b = 0; b < 3; ++b)
            std::cout << std::format(" {:+.4f}+/-{:.4f}", r.mc_cij(a, b),
                                     r.sigma_cij(a, b));
        std::cout << "\n";
    }
    std::cout << "\n Reconstructed B_+ (top polarization)\n";
    for (int a = 0; a < 3; ++a) {
        std::cout << std::format(
            " B+_{} : {:+.4f} +/- {:.4f}  (theory: {:+.4f})\n", ax[a],
            r.mc_bp(a), r.sigma_bp(a), wt.theory_bp(a));
    }

    std::cout << "\n Reconstructed B_- (anti-top polarization)\n";
    for (int a = 0; a < 3; ++a) {
        std::cout << std::format(
            " B-_{} : {:+.4f} +/- {:.4f}  (theory: {:+.4f})\n", ax[a],
            r.mc_bm(a), r.sigma_bm(a), wt.theory_bm(a));
    }

    std::cout << std::format(
        "\n Concurrence (MC)          : {:+.6f} +/- {:.6f}\n", r.mc_concurrence,
        r.sigma_concurrence);
    std::cout << std::format(" Concurrence (theory)      : {:+.6f}\n",
                             wt.theory_concurrence);
    std::cout << std::format(" significance(Concurrence) : {:.2f} sigma\n",
                             r.significance_concurrence);
    std::cout << std::format(
        "\n D=(Cnn-|Crr+Ckk|)/3 (MC)     : {:+.6f} +/- {:.6f}\n", r.mc_D,
        r.sigma_D);
    std::cout << std::format(" D=(Cnn-|Crr+Ckk|)/3 (theory) : {:+.6f}\n",
                             wt.theory_D);
    std::cout << std::format(" significance(D)              : {:.2f} sigma\n",
                             r.significance_D);

    std::cout << std::format("\n m12 (MC)         : {:.6f}\n", r.mc_m12);
    std::cout << std::format(" m12 - 1          : {:+.6f}\n", r.mc_m12 - 1.0);
    std::cout << std::format(" sigma_m12 (N_MC) : {:.6f}\n", r.sigma_m12);
    std::cout << std::format(" sig_m12 at N_MC : {:.2f} sigma\n",
                             r.significance_m12);

    // Phase 6: luminosity scan
    const std::vector<LumiScanPoint> lumi_scan =
        computeLumiScan(cfg, sigma_eff_fb, ev.n_accepted, r);

    // ------------------------------------------------------------------
    // Phase 7: fill and return MCResult
    // ------------------------------------------------------------------
    MCResult res;
    res.n_events_generated = ev.n_accepted;

    res.mc_bp = r.mc_bp;
    res.sigma_bp = r.sigma_bp;
    res.mc_bm = r.mc_bm;
    res.sigma_bm = r.sigma_bm;
    res.theory_bp = wt.theory_bp;
    res.theory_bm = wt.theory_bm;
    res.mc_cij = r.mc_cij;
    res.sigma_cij = r.sigma_cij;

    res.mc_D = r.mc_D;
    res.sigma_D = r.sigma_D;
    res.significance_D = r.significance_D;

    res.mc_concurrence = r.mc_concurrence;
    res.sigma_concurrence = r.sigma_concurrence;
    res.significance_concurrence = r.significance_concurrence;

    res.mc_m12 = r.mc_m12;
    res.sigma_m12 = r.sigma_m12;
    res.significance_m12 = r.significance_m12;

    res.theory_concurrence = wt.theory_concurrence;
    res.theory_D = wt.theory_D;
    res.theory_m12 = wt.theory_m12;

    res.total_xsec_fb = wt.total_weight;
    res.lumi_scan = std::move(lumi_scan);
    return res;
}
}  // namespace gagatt
