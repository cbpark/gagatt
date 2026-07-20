#ifndef SRC_MC_H
#define SRC_MC_H

#include <cstdint>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <Eigen/Dense>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include "constants.h"

namespace gagatt {
inline constexpr double X_DEFAULT = 4.8;  // machine parameter x

struct MCConfig {
    // Collider parameters
    double sqrt_s = 500.0;  // e+e- CM energy [GeV]
    double pe1 = +1.0;      // electron polarization Pe
    double pe2 = +1.0;      // positron polarization Pe~
    double x = X_DEFAULT;

    // Phase-space scan ranges
    double sqrts_min = 2.0 * MTOP + 1.0;  // just above threshold [GeV]
    double sqrts_max = -1.0;              // -1: use sqrt_s automatically
    double cos_th_min = -1.0;
    double cos_th_max = +1.0;

    // Grid for luminosity-weight z-cache
    int n_sqrts = 40;
    int n_cos = 20;

    // MC statistics
    long long n_events = 10'000'000LL;

    // RNG seed (0 = random_device)
    uint64_t seed = 42;

    // ----------------------------------------------------------------
    // Luminosity scan for significance plot
    // L_scan_min / L_scan_max  [ab^-1]   (0 disables the scan)
    // ----------------------------------------------------------------
    double L_scan_min_ab = 0.01;  // [ab^-1]
    double L_scan_max_ab = 0.2;   // [ab^-1]
    int n_L_points = 20;

    // Number of independent MC seeds to average over at each lumi point.
    // Set to 1 to reproduce the old single-seed behaviour (no spread).
    // Recommended: 10-20 for a smooth band; 1 is a fast approximation.
    int n_lumi_seeds = 10;
};

struct LumiScanPoint {
    double L_ab;  // luminosity [ab^-1]
    double significance_concurrence = 0.0;
    double significance_D = 0.0;    // (-1/3 - D) / sigma_D (0 if D >= -1/3)
    double significance_m12 = 0.0;  // (m12-1) / sigma_m12 (0 if m12 <= 1)

    // 1-sigma spread across the n_lumi_seeds realisations
    // (0 if n_lumi_seeds==1)
    double sigma_sig_concurrence = 0.0;
    double sigma_sig_D = 0.0;
    double sigma_sig_m12 = 0.0;
};

struct MCResult {
    long long n_events_generated = 0;

    // MC reconstructed quantities (at cfg.n_events)
    double mc_concurrence = 0.0;
    double mc_D = 0.0;
    double mc_m12 = 0.0;
    double sigma_concurrence = 0.0;
    double sigma_D = 0.0;
    double sigma_m12 = 0.0;

    // ----------------------------------------------------------------
    // Theory predictions (luminosity+phase-space weighted averages)
    // ----------------------------------------------------------------
    double theory_D = 0.0;
    double theory_concurrence = 0.0;
    double theory_m12 = 0.0;

    // Total cross section
    double total_xsec_fb = 0.0;  // [fb]

    // Luminosity scan (empty if L_scan_min_ab == 0)
    std::vector<LumiScanPoint> lumi_scan;
};

// Run the MC simulation and return aggregated results.
MCResult runMC(const MCConfig &cfg);

}  // namespace gagatt

#endif  // SRC_MC_H
