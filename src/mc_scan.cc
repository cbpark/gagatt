#include "mc_scan.h"
#include <format>
#include <iostream>
#include <random>
#include <vector>
#include "mc_helper.h"
#include "photon.h"
#include "spin_density.h"

namespace gagatt {
// -----------------------------------------------------------------------
// buildScanSliceWeightTable (internal helper)
//
// Builds a WeightTable sub-binned in cos(Theta) over [cos_th_min,
// cos_th_max] at a fixed sqrts_mid.
//
// The bin weight already absorbs d_sqrts_scan * d_cos_sub so that
// total_weight equals the integrated cross section [fb] over this slice.
// -----------------------------------------------------------------------
static WeightTable buildScanSliceWeightTable(
    double sqrts_mid, double d_sqrts_scan, const DScanConfig &scan_cfg,
    const LumiWeights &lw, double L_tot, const MCConfig &cfg) {
    const double d_cos_sub = (scan_cfg.cos_th_max - scan_cfg.cos_th_min) /
                             static_cast<double>(scan_cfg.n_cos_sub);

    WeightTable wt;
    wt.bin_weights.assign(scan_cfg.n_cos_sub, 0.0);
    wt.sdc_cache.reserve(scan_cfg.n_cos_sub);

    double tw_con = 0.0, tw_D = 0.0, tw_m12 = 0.0;
    for (int i = 0; i < scan_cfg.n_cos_sub; ++i) {
        const double cos_th = scan_cfg.cos_th_min + (i + 0.5) * d_cos_sub;
        const SDMatrixCoefficients sdc(sqrts_mid, cos_th, lw);

        const double rate = eventRate(sqrts_mid, sdc, L_tot, cfg.sqrt_s) *
                            d_sqrts_scan * d_cos_sub;

        wt.bin_weights[i] = std::max(0.0, rate);
        wt.sdc_cache.push_back(sdc);
        wt.total_weight += wt.bin_weights[i];

        tw_con += wt.bin_weights[i] * getConcurrence(sdc);
        tw_D += wt.bin_weights[i] * entanglementMarker(sdc);
        tw_m12 += wt.bin_weights[i] * horodeckiMeasure(sdc);
    }

    if (wt.total_weight > 0.0) {
        wt.theory_concurrence = tw_con / wt.total_weight;
        wt.theory_D = tw_D / wt.total_weight;
        wt.theory_m12 = tw_m12 / wt.total_weight;
    }

    return wt;
}

// -----------------------------------------------------------------------
// runDScanVsSqrtS
// -----------------------------------------------------------------------
std::vector<DScanBinResult> runDScanVsSqrtS(const MCConfig &cfg,
                                            const DScanConfig &scan_cfg,
                                            uint64_t seed) {
    const double sqrts_max = (cfg.sqrts_max > 0.0) ? cfg.sqrts_max : cfg.sqrt_s;
    const double sqrts_min = cfg.sqrts_min;

    const double pc1 = -cfg.pe1;  // PePc = -1
    const double pc2 = -cfg.pe2;

    const double d_sqrts_scan =
        (sqrts_max - sqrts_min) / static_cast<double>(scan_cfg.n_scan_bins);

    std::cout << std::format(
        "-- sqrt(s_hat) scan: {} bins in [{:.1f}, {:.1f}] GeV, "
        "cos(Theta) in [{:.2f}, {:.2f}], {} events/bin\n",
        scan_cfg.n_scan_bins, sqrts_min, sqrts_max, scan_cfg.cos_th_min,
        scan_cfg.cos_th_max, scan_cfg.n_events_per_bin);

    std::mt19937_64 rng(seed);

    std::vector<DScanBinResult> results;
    results.reserve(scan_cfg.n_scan_bins);

    for (int j = 0; j < scan_cfg.n_scan_bins; ++j) {
        const double sqrts_lo = sqrts_min + j * d_sqrts_scan;
        const double sqrts_hi = sqrts_lo + d_sqrts_scan;
        const double sqrts_mid = 0.5 * (sqrts_lo + sqrts_hi);

        DScanBinResult res;
        res.sqrts_lo = sqrts_lo;
        res.sqrts_hi = sqrts_hi;
        res.sqrts_mid = sqrts_mid;

        // Single lumi evaluation at the bin midpoint — valid when
        // d_sqrts_scan is small relative to the lumi variation scale.
        const double z = sqrts_mid / cfg.sqrt_s;
        const auto [lw, L_tot] =
            lumiWeightsAndTotal(z, cfg.x, cfg.pe1, pc1, cfg.pe2, pc2);

        const WeightTable wt = buildScanSliceWeightTable(
            sqrts_mid, d_sqrts_scan, scan_cfg, lw, L_tot, cfg);

        res.bin_xsec_fb = wt.total_weight;

        if (wt.total_weight <= 0.0) {
            // Below threshold or kinematically empty bin.
            results.push_back(res);
            continue;
        }

        res.theory_concurrence = wt.theory_concurrence;
        res.theory_D = wt.theory_D;
        res.theory_m12 = wt.theory_m12;

        // Dedicated MC event loop for this sqrt(s_hat) bin only.
        // verbose=false: suppress per-event progress prints across all bins.
        const EventLoopResult ev =
            runEventLoop(scan_cfg.n_events_per_bin, wt, rng, /*verbose=*/false);
        const ReconstructedMC r = reconstructFromMoments(ev);

        res.mc_concurrence = r.mc_concurrence;
        res.sigma_concurrence = r.sigma_concurrence;

        res.mc_D = r.mc_D;
        res.sigma_D = r.sigma_D;

        res.mc_m12 = r.mc_m12;
        res.sigma_m12 = r.sigma_m12;

        res.n_events = ev.n_accepted;

        results.push_back(res);
    }

    return results;
}
}  // namespace gagatt
