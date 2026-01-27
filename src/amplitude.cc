#include "amplitude.h"

#include <cmath>
#include <functional>
#include <numbers>

#include "constants.h"

#ifdef DEBUG
#include <iostream>
#endif

namespace gagatt {
constexpr double COUPLING_FACTOR = 8.0 * std::numbers::pi * ALPHA * QTOP2;

Amplitude offShellAmpApprox(double sqrt_s_hat, double cos_th, double m1,
                            double m2, Helicity lambda1, Helicity lambda2,
                            Helicity sigma1, Helicity sigma2) {
    // Input validation and threshold check
    const double threshold = m1 + m2;
    if (sqrt_s_hat < threshold || std::abs(cos_th) > 1.0) { return {0.0, 0.0}; }
    // r = 1 - beta^2
    const double s_hat = sqrt_s_hat * sqrt_s_hat;
    const double r = 2.0 * (m1 * m1 + m2 * m2) / s_hat;
    if (r > 1.0) { return {0.0, 0.0}; }

    // Kinematics
    const double beta2 = std::max(0.0, 1.0 - r);  // numerical safety
    const double beta = std::sqrt(beta2);
    const double sqrt_r = std::sqrt(r);  // sqrt(1 - beta2) = sqrt(r)
    const double cos_th2 = cos_th * cos_th;
    const double sin_th2 = std::max(0.0, 1.0 - cos_th2);  // numerical safety

    const double l1 = toDouble(lambda1);
    const double s1 = toDouble(sigma1);

    double amp = 0.0;

    if (lambda1 == lambda2) {
        if (sigma1 == sigma2) { amp = sqrt_r * (beta * s1 + l1); }
    } else {  // lambda1 == - lambda2
        if (sigma1 == sigma2) {
            amp = -beta * sqrt_r * sin_th2 * s1;
        } else {  // lambda1 == - lambda2 && sigma1 == -sigma2
            const double sin_th = std::sqrt(sin_th2);
            amp = -beta * sin_th * (s1 * l1 + cos_th);
        }
    }

    // Early exit for forbidden helicity combinations
    if (amp == 0.0) { return {0.0, 0.0}; }

    // : sin^2 + r*cos^2 is more stable than 1 - beta^2*cos^2
    const double denom = sin_th2 + r * cos_th2;
    if (denom < 1e-18) { return {0.0, 0.0}; }

    return {(COUPLING_FACTOR / denom) * amp, 0.0};
}

double c1OnShell(double sqrt_s_hat, double cos_th) {
    return lam1lam2Sum([=](Helicity l1, Helicity l2) {
        using H = Helicity;
        return onShellHelAmp2(sqrt_s_hat, cos_th, l1, l2, H::PLUS, H::PLUS) +
               onShellHelAmp2(sqrt_s_hat, cos_th, l1, l2, H::MINUS, H::MINUS);
    });
}

double c2OnShell(double sqrt_s_hat, double cos_th) {
    return lam1lam2Sum([=](Helicity l1, Helicity l2) {
        using H = Helicity;
        return onShellHelAmp2(sqrt_s_hat, cos_th, l1, l2, H::PLUS, H::PLUS) -
               onShellHelAmp2(sqrt_s_hat, cos_th, l1, l2, H::MINUS, H::MINUS);
    });
}

double c3OnShell(double sqrt_s_hat, double cos_th) {
    return lam1lam2Sum([=](Helicity l1, Helicity l2) {
        using H = Helicity;
        return onShellHelAmp2(sqrt_s_hat, cos_th, l1, l2, H::MINUS, H::PLUS) +
               onShellHelAmp2(sqrt_s_hat, cos_th, l1, l2, H::PLUS, H::MINUS);
    });
}

double c1OffShellApprox(double sqrt_s_hat, double cos_th, double m1,
                        double m2) {
    return lam1lam2Sum([=](Helicity l1, Helicity l2) {
        using H = Helicity;
        return offShellHelAmp2Approx(sqrt_s_hat, cos_th, m1, m2, l1, l2,
                                     H::PLUS, H::PLUS) +
               offShellHelAmp2Approx(sqrt_s_hat, cos_th, m1, m2, l1, l2,
                                     H::MINUS, H::MINUS);
    });
}
}  // namespace gagatt
