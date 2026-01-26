#ifndef SRC_AMPLITUDE_H
#define SRC_AMPLITUDE_H

#include <complex>
#include "constants.h"

namespace gagatt {
using Amplitude = std::complex<double>;

enum class Polarization { PLUS, MINUS };

inline constexpr double toDouble(Polarization pol) {
    return (pol == Polarization::PLUS) ? 1.0 : -1.0;
}

inline constexpr Polarization operator-(Polarization pol) {
    return (pol == Polarization::PLUS) ? Polarization::MINUS
                                       : Polarization::PLUS;
}

Amplitude offShellAmplitudeApprox(double s_hat, double cos_th, double m1,
                                  double m2, Polarization lambda1,
                                  Polarization lambda2, Polarization sigma1,
                                  Polarization sigma2);

inline Amplitude onShellAmplitude(double s_hat, double cos_th,
                                  Polarization lambda1, Polarization lambda2,
                                  Polarization sigma1, Polarization sigma2) {
    return offShellAmplitudeApprox(s_hat, cos_th, MTOP, MTOP, lambda1, lambda2,
                                   sigma1, sigma2);
}

inline double onShellHelAmp2(double s_hat, double cos_th, Polarization lambda1,
                             Polarization lambda2, Polarization sigma1,
                             Polarization sigma2) {
    return std::norm(
        onShellAmplitude(s_hat, cos_th, lambda1, lambda2, sigma1, sigma2));
}

template <typename F>
constexpr auto lam1lam2Sum(F &&f)
    -> decltype(f(Polarization::PLUS, Polarization::PLUS)) {
    using P = Polarization;
    return 0.25 * (f(P::PLUS, P::PLUS) + f(P::PLUS, P::MINUS) +
                   f(P::MINUS, P::PLUS) + f(P::MINUS, P::MINUS));
}

double c1OnShell(double s_hat, double cos_th);
double c2OnShell(double s_hat, double cos_th);
double c3OnShell(double s_hat, double cos_th);

double onShellAmp2Sum(double s_hat, double cos_th) {
    return c1OnShell(s_hat, cos_th) + c3OnShell(s_hat, cos_th);
}
}  // namespace gagatt

#endif  // SRC_AMPLITUDE_H
