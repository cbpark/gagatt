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
    double L_ee_fb = 1000.0;  // integrated luminosity [fb^-1]

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
};

struct MCResult {
    long long n_events_generated = 0;

    // ----------------------------------------------------------------
    // Reconstructed C_ij from <q+_i q-_j> = -1/9 C_ij
    // ----------------------------------------------------------------
    Eigen::Matrix3d mc_cij = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d sigma_cij = Eigen::Matrix3d::Zero();

    // Tr[C] and the entanglement marker D = Tr[C]/3
    double mc_tr_c = 0.0;
    double sigma_tr_c = 0.0;

    // Significance of D < -1/3 vs null (D = 0): |D| / sigma_D
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
    double theory_tr_c = 0.0;
    double theory_concurrence = 0.0;
    double theory_negativity = 0.0;
    double theory_m12 = 0.0;

    // Total cross section
    double total_xsec_fb = 0.0;  // [fb]
};

// Run the MC simulation and return aggregated results.
// Prints progress to stdout.
MCResult runMC(const MCConfig &cfg);

}  // namespace gagatt

#endif  // SRC_MC_H
