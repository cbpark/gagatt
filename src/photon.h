#ifndef SRC_PHOTON_H
#define SRC_PHOTON_H

#include <array>
#include <optional>
#include "helicity.h"

namespace gagatt {
double fLumi(double x, double y, double pe, double pc);

// z = sqrt(tau)
double photonLuminosity(double z, double x, double pe1, double pc1, double pe2,
                        double pc2, std::optional<Helicity> l1,
                        std::optional<Helicity> l2);

// For unpoloarized initial beams (pe1 = pc1 = pe2 = pc2 = 0)
// In this case, only contribution from <00> remains.
// Note that, 0.25 * sum_{l1, l2} <00> = <00> with l1 = l2 = 0.
inline double photonLuminosityUnpol(double z, double x) {
    return photonLuminosity(z, x, 0.0, 0.0, 0.0, 0.0, {}, {});
}

// Normalised luminosity weights for the four photon helicity combinations.
// Order: {++, +-, -+, --}.
struct LumiWeights {
    std::array<double, 4> w;  // democratic default
};

// Compute luminosity weights from polarization correlation integrals.
// Performs 4 GSL integrations (c00, c20, c02, c22) and returns the
// normalised weights for (++, +-, -+, --).
LumiWeights lumiWeights(double z, double x, double pe1, double pc1, double pe2,
                        double pc2);

// Luminosity weights plus the un-normalised total L_tot = sum_{l1,l2} L^{l1 l2}
// L_tot has units of [1/GeV^2] (same as each individual L^{l1 l2}).
// The MC needs this to compute the absolute event rate.
struct LumiWeightsAndTotal {
    LumiWeights lw;
    double L_tot{0.0};  // sum_{l1,l2} L^{l1 l2}(z)  [GeV^-2]
};

LumiWeightsAndTotal lumiWeightsAndTotal(double z, double x, double pe1,
                                        double pc1, double pe2, double pc2);
}  // namespace gagatt

#endif  // SRC_PHOTON_H
