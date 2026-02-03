#ifndef SRC_AMPLITUDE_H
#define SRC_AMPLITUDE_H

#include <array>
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

inline constexpr std::array<Helicity, 2> all_helicities = {Helicity::PLUS,
                                                           Helicity::MINUS};

template <typename F>
constexpr auto averageHelicities(F &&f) {
    auto sum = f(Helicity::PLUS, Helicity::PLUS) * 0.0;  // deduce type
    for (auto h1 : all_helicities) {
        for (auto h2 : all_helicities) { sum += f(h1, h2); }
    }
    return sum * 0.25;
}

double c1OnShell(double sqrt_s_hat, double cos_th);
double c2OnShell(double sqrt_s_hat, double cos_th);
double c3OnShell(double sqrt_s_hat, double cos_th);
double c4OnShell(double sqrt_s_hat, double cos_th);

inline double onShellAmp2Sum(double sqrt_s_hat, double cos_th) {
    return c1OnShell(sqrt_s_hat, cos_th) + c3OnShell(sqrt_s_hat, cos_th);
}

// double c1OffShellApprox(double sqrt_s_hat, double cos_th, double m1, double
// m2); double c1TildeOffShellApprox(double sqrt_s_hat, double cos_th);
}  // namespace gagatt

#endif  // SRC_AMPLITUDE_H
