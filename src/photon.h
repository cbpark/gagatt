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
}  // namespace gagatt

#endif  // SRC_PHOTON_H
