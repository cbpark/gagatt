#ifndef SRC_AMPLITUDE_H
#define SRC_AMPLITUDE_H

#include <complex>
#include "constants.h"

namespace gagatt {
using Amplitude = std::complex<double>;

enum class Helicity { PLUS, MINUS };

inline constexpr double toDouble(Helicity pol) {
    return (pol == Helicity::PLUS) ? 1.0 : -1.0;
}

inline constexpr Helicity operator-(Helicity pol) {
    return (pol == Helicity::PLUS) ? Helicity::MINUS : Helicity::PLUS;
}

Amplitude offShellAmplitudeApprox(double s_hat, double cos_th, double m1,
                                  double m2, Helicity lambda1, Helicity lambda2,
                                  Helicity sigma1, Helicity sigma2);

inline Amplitude onShellAmplitude(double s_hat, double cos_th, Helicity lambda1,
                                  Helicity lambda2, Helicity sigma1,
                                  Helicity sigma2) {
    return offShellAmplitudeApprox(s_hat, cos_th, MTOP, MTOP, lambda1, lambda2,
                                   sigma1, sigma2);
}

inline double onShellHelAmp2(double s_hat, double cos_th, Helicity lambda1,
                             Helicity lambda2, Helicity sigma1,
                             Helicity sigma2) {
    return std::norm(
        onShellAmplitude(s_hat, cos_th, lambda1, lambda2, sigma1, sigma2));
}

template <typename F>
constexpr auto lam1lam2Sum(F &&f)
    -> decltype(f(Helicity::PLUS, Helicity::PLUS)) {
    using H = Helicity;
    return 0.25 * (f(H::PLUS, H::PLUS) + f(H::PLUS, H::MINUS) +
                   f(H::MINUS, H::PLUS) + f(H::MINUS, H::MINUS));
}

double c1OnShell(double s_hat, double cos_th);
double c2OnShell(double s_hat, double cos_th);
double c3OnShell(double s_hat, double cos_th);

double onShellAmp2Sum(double s_hat, double cos_th) {
    return c1OnShell(s_hat, cos_th) + c3OnShell(s_hat, cos_th);
}
}  // namespace gagatt

#endif  // SRC_AMPLITUDE_H
