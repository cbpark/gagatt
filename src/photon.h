#ifndef SRC_PHOTON_H
#define SRC_PHOTON_H

namespace gagatt {
double fLumi(double x, double y, double pe, double pc);

// sigma_c * 2 pi alpha^2 / (x m_e^2)
double sigmaC(double x, double pe, double pc);

// struct PhotonPolarization {
//     double c0, c2;
// };

// PhotonPolarization computePhotonC(double x, double y);
}  // namespace gagatt

#endif  // SRC_PHOTON_H
