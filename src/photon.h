#ifndef SRC_PHOTON_H
#define SRC_PHOTON_H

#include <array>
#include <optional>
#include "helicity.h"

namespace gagatt {
double fLumi(double x, double y, double pe, double pc);

// Normalised luminosity weights for the four photon helicity combinations.
// Order: {++, +-, -+, --}.
struct LumiWeights {
    std::array<double, 4> w;  // normalised: sum = 1
};

// Luminosity weights plus the un-normalised total L_tot = sum_{l1,l2} L^{l1 l2}
// L_tot has units of [1/GeV^2].  This is the primary computation; all other
// lumi functions are implemented in terms of this one.
struct LumiWeightsAndTotal {
    LumiWeights lw;
    double L_tot{0.0};  // sum_{l1,l2} L^{l1 l2}(z)  [GeV^-2]
};

// Primary function: computes {LumiWeights, L_tot} in a single pass.
// Performs 4 GSL integrations (c00, c20, c02, c22).
// z = sqrt(tau) = sqrt_s_hat / sqrt_s
LumiWeightsAndTotal lumiWeightsAndTotal(double z, double x, double pe1,
                                        double pc1, double pe2, double pc2);

// Convenience wrapper: returns only the normalised weights.
inline LumiWeights lumiWeights(double z, double x, double pe1, double pc1,
                               double pe2, double pc2) {
    return lumiWeightsAndTotal(z, x, pe1, pc1, pe2, pc2).lw;
}

// Full polarised photon luminosity for a specific helicity pair (l1, l2).
// Pass {} for an unpolarised (averaged) helicity.
// z = sqrt(tau)
double photonLuminosity(double z, double x, double pe1, double pc1, double pe2,
                        double pc2, std::optional<Helicity> l1,
                        std::optional<Helicity> l2);

// Unpolarised photon luminosity (pe1 = pc1 = pe2 = pc2 = 0).
// Note: 0.25 * sum_{l1,l2} L^{l1 l2} = L with l1 = l2 = 0.
inline double photonLuminosityUnpol(double z, double x) {
    return photonLuminosity(z, x, 0.0, 0.0, 0.0, 0.0, {}, {});
}
}  // namespace gagatt

#endif  // SRC_PHOTON_H
