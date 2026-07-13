#ifndef SRC_MC_H
#define SRC_MC_H

#include <cstdint>
#include "constants.h"

namespace gagatt {
inline constexpr double X_DEFAULT = 4.8;  // machine parameter x
// analyzing power for charged leptons
inline constexpr double KAPPA_LEPTON = 1.0;

struct MCConfig {
    // Collider
    double sqrt_s = 500.0;  // e+e- CM energy [GeV]
    double pe1 = +1.0;      // electron polarization
    double pe2 = +1.0;      // positron polarization
    double x = X_DEFAULT;
    double L_ee_fb = 1000.0;  // integrated luminosity [fb^{-1}]

    // Scan ranges
    double sqrts_min = 2.0 * MTOP + 1.0;  // just above threshold [GeV]
    double sqrts_max = -1.0;              // -1: set to sqrt_s automatically
    double cos_th_min = -1.0;
    double cos_th_max = +1.0;

    // Grid for the z-cache (lumi weights precomputation)
    int n_sqrts = 2000;
    int n_cos = 200;

    // MC statistics
    long long n_events = 1'000'000LL;

    // RNG seed (0 = random device)
    uint64_t seed = 42;
};

struct MCResult {
    long long n_events_generated = 0;  // total accepted events
    double mean_cos_phi = 0;           // <cos phi>
    double sigma_cos_phi = 0;          // statistical uncertainty
    double significance = 0;  // |<cos phi>| / sigma  (null: <cos phi>=0)

    // Theoretical prediction from the density matrix
    double theory_cos_phi = 0;      // -Tr[C]/9
    double theory_negativity = 0;   // N[rho] integrated over phase space
    double theory_concurrence = 0;  // C[rho] integrated over phase space

    // Cross section
    double total_xsec_fb = 0;  // total cross section [fb]
};

// Run the MC and return aggregated results.
// Prints progress to stdout.
MCResult runMC(const MCConfig &cfg);
}  // namespace gagatt

#endif  // SRC_MC_H
