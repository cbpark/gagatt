#ifndef SRC_PHOTON_H
#define SRC_PHOTON_H

#include <optional>
#include "amplitude.h"

namespace gagatt {
double fLumi(double x, double y, double pe, double pc);

// z = sqrt(tau)
double photonLuminosity(double z, double x, double pe1, double pc1, double pe2,
                        double pc2, std::optional<Helicity> l1,
                        std::optional<Helicity> l2);

inline double photonLuminosityUnpol(double z, double x, double pe1, double pc1,
                             double pe2, double pc2) {
    return photonLuminosity(z, x, pe1, pc1, pe2, pc2, {}, {});
}
}  // namespace gagatt

#endif  // SRC_PHOTON_H
