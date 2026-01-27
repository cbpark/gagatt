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

Amplitude offShellAmpApprox(double sqrt_s_hat, double cos_th, double m1,
                            double m2, Helicity lambda1, Helicity lambda2,
                            Helicity sigma1, Helicity sigma2);

inline double offShellHelAmp2Approx(double sqrt_s_hat, double cos_th, double m1,
                                    double m2, Helicity lambda1,
                                    Helicity lambda2, Helicity sigma1,
                                    Helicity sigma2) {
    return std::norm(offShellAmpApprox(sqrt_s_hat, cos_th, m1, m2, lambda1,
                                       lambda2, sigma1, sigma2));
}

inline Amplitude onShellAmp(double sqrt_s_hat, double cos_th, Helicity lambda1,
                            Helicity lambda2, Helicity sigma1,
                            Helicity sigma2) {
    return offShellAmpApprox(sqrt_s_hat, cos_th, MTOP, MTOP, lambda1, lambda2,
                             sigma1, sigma2);
}

inline double onShellHelAmp2(double sqrt_s_hat, double cos_th, Helicity lambda1,
                             Helicity lambda2, Helicity sigma1,
                             Helicity sigma2) {
    return offShellHelAmp2Approx(sqrt_s_hat, cos_th, MTOP, MTOP, lambda1,
                                 lambda2, sigma1, sigma2);
}

template <typename M2>
constexpr auto lam1lam2Sum(M2 &&f)
    -> decltype(f(Helicity::PLUS, Helicity::PLUS)) {
    using H = Helicity;
    return 0.25 * (f(H::PLUS, H::PLUS) + f(H::PLUS, H::MINUS) +
                   f(H::MINUS, H::PLUS) + f(H::MINUS, H::MINUS));
}

double c1OnShell(double sqrt_s_hat, double cos_th);
double c2OnShell(double sqrt_s_hat, double cos_th);
double c3OnShell(double sqrt_s_hat, double cos_th);

inline double onShellAmp2Sum(double sqrt_s_hat, double cos_th) {
    return c1OnShell(sqrt_s_hat, cos_th) + c3OnShell(sqrt_s_hat, cos_th);
}

double c1OffShellApprox(double sqrt_s_hat, double cos_th, double m1, double m2);
}  // namespace gagatt

#endif  // SRC_AMPLITUDE_H
