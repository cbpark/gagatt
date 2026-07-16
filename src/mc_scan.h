#ifndef SRC_MC_SCAN_H
#define SRC_MC_SCAN_H

#include <cstdint>
#include <vector>
#include "mc.h"
#include "mc_helper.h"

namespace gagatt {
// -----------------------------------------------------------------------
// DScanConfig
//
// Configuration for a 1-D scan of observable vs sqrt(s_hat), with a fixed
// cut on cos(Theta) applied within each sqrt(s_hat) bin (cos(Theta) is
// integrated/sampled only over [cos_th_min, cos_th_max], not scanned).
// -----------------------------------------------------------------------
struct DScanConfig {
    int n_scan_bins = 50;                 // number of sqrt(s_hat) scan bins
    int n_cos_sub = 50;                   // cos(Theta) sub-bins within each
                                          //   sqrt(s_hat) bin
    long long n_events_per_bin = 100000;  // MC events per sqrt(s_hat) bin
    double cos_th_min = -1.0;             // cos(Theta) cut window
    double cos_th_max = 1.0;
};

// -----------------------------------------------------------------------
// DScanBinResult
//
// Observable quantities in one sqrt(s_hat) bin,
// integrated/sampled over the fixed cos(Theta) cut window.
// -----------------------------------------------------------------------
struct DScanBinResult {
    double sqrts_lo = 0.0;
    double sqrts_hi = 0.0;
    double sqrts_mid = 0.0;

    double theory_concurrence = 0.0;
    double theory_D = 0.0;

    double mc_concurrence = 0.0;
    double sigma_concurrence = 0.0;
    double significance_concurrence = 0.0;

    double mc_D = 0.0;
    double sigma_D = 0.0;
    double significance_D = 0.0;

    double bin_xsec_fb = 0.0;  // integrated xsec over cos_th cut window
    long long n_events = 0;
};

// -----------------------------------------------------------------------
// runDScanVsSqrtS
//
// Scans D and concurrence vs sqrt(s_hat) with a fixed cos(Theta) cut
// window defined in scan_cfg.  For each sqrt(s_hat) bin:
//   - evaluates the lumi-weighted theory values at the bin midpoint
//   - runs a dedicated MC event loop (n_events_per_bin events)
//   - collects DScanBinResult (theory/MC values + significance)
//
// Sub-binning in cos(Theta) uses scan_cfg.n_cos_sub bins over
// [cos_th_min, cos_th_max].  A single lumi evaluation at the bin midpoint
// is used (valid when d_sqrts_scan is small enough that the lumi is
// smooth over one bin width).
// -----------------------------------------------------------------------
std::vector<DScanBinResult> runDScanVsSqrtS(const MCConfig &cfg,
                                            const DScanConfig &scan_cfg,
                                            uint64_t seed);
}  // namespace gagatt

#endif  // SRC_MC_SCAN_H
