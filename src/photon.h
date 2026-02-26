#ifndef SRC_PHOTON_H
#define SRC_PHOTON_H

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

struct LumiWeights {
    double wpp, wmm, wpm, wmp;
};

// Compute the four luminosity weights w_{l1l2} = L_{l1l2} / sum L
// at a given z = sqrt(tau).
// Returns {0,0,0,0} if sum vanishes.
LumiWeights lumiWeights(double z, double x, double pe1, double pc1, double pe2,
                        double pc2);
}  // namespace gagatt

#endif  // SRC_PHOTON_H
