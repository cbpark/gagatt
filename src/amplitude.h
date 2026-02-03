#ifndef SRC_AMPLITUDE_H
#define SRC_AMPLITUDE_H

#include <array>
#include <complex>
#include <functional>
#include <type_traits>
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

struct PolarizationCoefficients {
    double c1, c2, c3, c4;
};

inline constexpr std::array<Helicity, 2> all_helicities = {Helicity::PLUS,
                                                           Helicity::MINUS};

template <typename F>
constexpr auto averageHelicities(F &&f) {
    using ResultType = std::invoke_result_t<F, Helicity, Helicity>;
    ResultType total{};

    for (auto l1 : all_helicities) {
        for (auto l2 : all_helicities) {
            auto res = f(l1, l2);
            if constexpr (std::is_arithmetic_v<ResultType>) {
                total += res;
            } else {
                total.c1 += res.c1;
                total.c2 += res.c2;
                total.c3 += res.c3;
                total.c4 += res.c4;
            }
        }
    }

    if constexpr (std::is_arithmetic_v<ResultType>) {
        return total * 0.25;
    } else {
        return PolarizationCoefficients{total.c1 * 0.25, total.c2 * 0.25,
                                        total.c3 * 0.25, total.c4 * 0.25};
    }
}

PolarizationCoefficients computeCoeffsOnShell(double sqrt_s_hat, double cos_th);
}  // namespace gagatt

#endif  // SRC_AMPLITUDE_H
