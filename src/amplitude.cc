#include "amplitude.h"

#include <cmath>
#include <functional>
#include <numbers>

#include "constants.h"

namespace gagatt {
constexpr double COUPLING_FACTOR = 8.0 * std::numbers::pi * ALPHA * QTOP2;

Amplitude offShellAmplitudeApprox(const double s_hat, const double cos_th,
                                  const double m1, const double m2,
                                  const Polarization lambda1,
                                  const Polarization lambda2,
                                  const Polarization sigma1,
                                  const Polarization sigma2) {
    // Input validation and threshold check
    if (cos_th > 1.0 || cos_th < -1.0) { return 0.0; }

    const double ratio = 2.0 * (m1 * m1 + m2 * m2) / s_hat;
    if (ratio > 1.0) { return 0.0; }

    // Kinematics
    const double beta2 = 1.0 - ratio;
    const double beta = std::sqrt(beta2);
    // Simplified: sqrt(1 - beta2) = sqrt(ratio)
    const double sqrt_ratio = std::sqrt(ratio);
    const double cos_th2 = cos_th * cos_th;
    const double sin_th2 = std::max(0.0, 1.0 - cos_th2);  // Numerical safety

    const double l1 = toDouble(lambda1);
    const double s1 = toDouble(sigma1);

    double amp = 0.0;

    if (lambda1 == lambda2) {
        if (sigma1 == sigma2) { amp = sqrt_ratio * (beta * s1 + l1); }
    } else {  // lambda1 == - lambda2
        if (sigma1 == sigma2) {
            amp = -beta * sqrt_ratio * sin_th2;
        } else {  // lambda1 == - lambda2 && sigma1 == -sigma2
            const double sin_th = std::sqrt(sin_th2);
            amp = -beta * sin_th * (s1 * l1 + cos_th);
        }
    }

    // Early exit if zero
    if (amp == 0.0) { return 0.0; }

    // : sin^2 + ratio*cos^2 is more stable than 1 - beta^2*cos^2
    const double denom = sin_th2 + ratio * cos_th2;

    return {(COUPLING_FACTOR / denom) * amp, 0.0};
}

template <typename F>
constexpr auto lam1lam2Sum(F &&f)
    -> decltype(f(Polarization::PLUS, Polarization::PLUS)) {
    using P = Polarization;

    return 0.25 * (f(P::PLUS, P::PLUS) + f(P::PLUS, P::MINUS) +
                   f(P::MINUS, P::PLUS) + f(P::MINUS, P::MINUS));
}

inline double onShellAmp2(double s, double c, Polarization l1, Polarization l2,
                          Polarization s1, Polarization s2) {
    return std::norm(onShellAmplitude(s, c, l1, l2, s1, s2));
}

double c1OnShell(const double s_hat, const double cos_th) {
    return lam1lam2Sum([=](Polarization l1, Polarization l2) {
        using P = Polarization;
        return onShellAmp2(s_hat, cos_th, l1, l2, P::PLUS, P::PLUS) +
               onShellAmp2(s_hat, cos_th, l1, l2, P::MINUS, P::MINUS);
    });
}

double c2OnShell(const double s_hat, const double cos_th) {
    return lam1lam2Sum([=](Polarization l1, Polarization l2) {
        using P = Polarization;
        return onShellAmp2(s_hat, cos_th, l1, l2, P::PLUS, P::PLUS) -
               onShellAmp2(s_hat, cos_th, l1, l2, P::MINUS, P::MINUS);
    });
}

double c3OnShell(const double s_hat, const double cos_th) {
    return lam1lam2Sum([=](Polarization l1, Polarization l2) {
        using P = Polarization;
        return onShellAmp2(s_hat, cos_th, l1, l2, P::MINUS, P::PLUS) +
               onShellAmp2(s_hat, cos_th, l1, l2, P::PLUS, P::MINUS);
    });
}
}  // namespace gagatt
