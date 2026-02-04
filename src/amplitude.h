#ifndef SRC_AMPLITUDE_H
#define SRC_AMPLITUDE_H

#include <complex>
#include <type_traits>
#include "constants.h"
#include "helicity.h"

namespace gagatt {
using Amplitude = std::complex<double>;

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

inline Amplitude onShellHelAmp(double sqrt_s_hat, double cos_th,
                               Helicity lambda1, Helicity lambda2,
                               Helicity sigma1, Helicity sigma2) {
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
    double c1{0}, c2{0}, c3{0}, c4{0}, c5{0}, c6{0}, c7{0}, c8{0}, c9{0},
        c10{0}, c11{0}, c12{0}, c13{0}, c14{0}, c15{0}, c16{0};

    PolarizationCoefficients &operator+=(
        const PolarizationCoefficients &other) {
        c1 += other.c1;
        c2 += other.c2;
        c3 += other.c3;
        c4 += other.c4;
        c5 += other.c5;
        c6 += other.c6;
        c7 += other.c7;
        c8 += other.c8;
        c9 += other.c9;
        c10 += other.c10;
        c11 += other.c11;
        c12 += other.c12;
        c13 += other.c13;
        c14 += other.c14;
        c15 += other.c15;
        c16 += other.c16;
        return *this;
    }

    PolarizationCoefficients operator*(double factor) const {
        return {c1 * factor,  c2 * factor,  c3 * factor,  c4 * factor,
                c5 * factor,  c6 * factor,  c7 * factor,  c8 * factor,
                c9 * factor,  c10 * factor, c11 * factor, c12 * factor,
                c13 * factor, c14 * factor, c15 * factor, c16 * factor};
    }
};

template <typename F>
auto averageHelicities(F &&func) {
    using ResultType = std::invoke_result_t<F, Helicity, Helicity>;
    ResultType total{};

    for (auto l1 : {Helicity::PLUS, Helicity::MINUS}) {
        for (auto l2 : {Helicity::PLUS, Helicity::MINUS}) {
            total += func(l1, l2);
        }
    }

    return total * 0.25;
}

PolarizationCoefficients computeCoeffsOnShell(double sqrt_s_hat, double cos_th);
}  // namespace gagatt

#endif  // SRC_AMPLITUDE_H
