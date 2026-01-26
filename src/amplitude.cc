#include "amplitude.h"

#include <cmath>
#include <functional>
#include <numbers>

#include "constants.h"

namespace gagatt {
constexpr double COUPLING_FACTOR = 8.0 * std::numbers::pi * ALPHA * QTOP2;

Amplitude offShellAmplitudeApprox(double s_hat, double cos_th, double m1,
                                  double m2, Polarization lambda1,
                                  Polarization lambda2, Polarization sigma1,
                                  Polarization sigma2) {
    // Input validation and threshold check
    if (cos_th > 1.0 || cos_th < -1.0) { return 0.0; }

    const double ratio = 2.0 * (m1 * m1 + m2 * m2) / s_hat;
    if (ratio > 1.0) { return 0.0; }

    // Kinematics
    const double beta2 = 1.0 - ratio;
    const double beta = std::sqrt(beta2);
    const double sqrt_ratio =
        std::sqrt(ratio);  // sqrt(1 - beta2) = sqrt(ratio)
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

double c1OnShell(double s_hat, double cos_th) {
    return lam1lam2Sum([=](Polarization l1, Polarization l2) {
        using P = Polarization;
        return onShellHelAmp2(s_hat, cos_th, l1, l2, P::PLUS, P::PLUS) +
               onShellHelAmp2(s_hat, cos_th, l1, l2, P::MINUS, P::MINUS);
    });
}

double c2OnShell(double s_hat, double cos_th) {
    return lam1lam2Sum([=](Polarization l1, Polarization l2) {
        using P = Polarization;
        return onShellHelAmp2(s_hat, cos_th, l1, l2, P::PLUS, P::PLUS) -
               onShellHelAmp2(s_hat, cos_th, l1, l2, P::MINUS, P::MINUS);
    });
}

double c3OnShell(double s_hat, double cos_th) {
    return lam1lam2Sum([=](Polarization l1, Polarization l2) {
        using P = Polarization;
        return onShellHelAmp2(s_hat, cos_th, l1, l2, P::MINUS, P::PLUS) +
               onShellHelAmp2(s_hat, cos_th, l1, l2, P::PLUS, P::MINUS);
    });
}
}  // namespace gagatt
