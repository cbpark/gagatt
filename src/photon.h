#ifndef SRC_PHOTON_H
#define SRC_PHOTON_H

namespace gagatt {
double fLumi(double x, double y, double pe, double pc);

// z = sqrt(tau)
double photonLuminosityUnpol(double z, double x, double pe1, double pc1,
                             double pe2, double pc2);
}  // namespace gagatt

#endif  // SRC_PHOTON_H
