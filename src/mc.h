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
inline constexpr double X_DEFAULT = 4.8;     // machine parameter x
inline constexpr double KAPPA_LEPTON = 1.0;  // analyzing power for l+/-

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
    int n_sqrts = 2000;
    int n_cos = 200;

    // MC statistics
    long long n_events = 1'000'000LL;

    // RNG seed (0 = random_device)
    uint64_t seed = 42;

    // ----------------------------------------------------------------
    // Luminosity scan for significance plot
    // L_scan_min / L_scan_max  [ab^-1]   (0 disables the scan)
    // ----------------------------------------------------------------
    double L_scan_min_ab = 0.01;  // [ab^-1]
    double L_scan_max_ab = 2.0;   // [ab^-1]
    int n_L_points = 200;
};

struct LumiScanPoint {
    double L_ab;               // luminosity [ab^-1]
    double significance_D;     // (|D| - 1/3) / sigma_D (0 if D > -1/3)
    double significance_bell;  // (m12-1) / sigma_m12  (0 if m12 <= 1)
};

struct MCResult {
    long long n_events_generated = 0;

    // ----------------------------------------------------------------
    // Reconstructed C_ij from <q+_i q-_j> = -1/9 C_ij
    // ----------------------------------------------------------------
    Eigen::Matrix3d mc_cij = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d sigma_cij = Eigen::Matrix3d::Zero();
    Eigen::Vector3d mc_bp = Eigen::Vector3d::Zero();
    Eigen::Vector3d sigma_bp = Eigen::Vector3d::Zero();
    Eigen::Vector3d mc_bm = Eigen::Vector3d::Zero();
    Eigen::Vector3d sigma_bm = Eigen::Vector3d::Zero();

    // Tr[C]. note that D =/= Tr[C] / 3
    double mc_tr_c = 0.0;
    double sigma_tr_c = 0.0;

    // Entanglement marker: D = (C_nn - |C_rr + C_kk|) / 3
    double mc_D = 0.0;
    double sigma_D = 0.0;

    // Significance of entanglement: (-D - 1/3) / sigma_D when D < -1/3
    double significance_D = 0.0;

    // Quantum-information quantities derived from reconstructed C_ij
    // (density matrix built with B+ = B- = 0, valid at LO in QED)
    double mc_concurrence = 0.0;
    double mc_negativity = 0.0;
    double mc_m12 = 0.0;  // Horodecki m12 = m1 + m2

    // Significance of Bell inequality violation: (m12 - 1) / sigma_m12
    double significance_bell = 0.0;

    // ----------------------------------------------------------------
    // Theory predictions (luminosity+phase-space weighted averages)
    // ----------------------------------------------------------------
    double theory_tr_c = 0.0;  // Tr[C]
    double theory_D = 0.0;     // (C_nn - |C_rr + C_kk|) / 3
    double theory_concurrence = 0.0;
    double theory_negativity = 0.0;
    double theory_m12 = 0.0;

    // Total cross section
    double total_xsec_fb = 0.0;  // [fb]

    // Luminosity scan (empty if L_scan_min_ab == 0)
    std::vector<LumiScanPoint> lumi_scan;
};

// Run the MC simulation and return aggregated results.
// Prints progress to stdout.
MCResult runMC(const MCConfig &cfg);

}  // namespace gagatt

#endif  // SRC_MC_H
