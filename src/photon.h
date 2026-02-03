#ifndef SRC_PHOTON_H
#define SRC_PHOTON_H

#include "amplitude.h"

namespace gagatt {
// sigma_c * 2 pi alpha^2 / (x m_e^2)
double sigmaC(double x, Helicity pe, Helicity pc);
}  // namespace gagatt

#endif  // SRC_PHOTON_H
