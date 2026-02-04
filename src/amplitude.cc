#include "amplitude.h"

#include <cmath>
#include <complex>
#include <functional>

#include "constants.h"
#include "helicity.h"

#ifdef DEBUG
#include <iostream>
#endif

namespace gagatt {
Amplitude offShellAmpApprox(double sqrt_s_hat, double cos_th, double m1,
                            double m2, Helicity lambda1, Helicity lambda2,
                            Helicity sigma1, Helicity sigma2) {
    // Input validation and threshold check
    const double threshold = m1 + m2;
    if (sqrt_s_hat < threshold || std::abs(cos_th) > 1.0) {
#ifdef DEBUG
        std::cerr << "offShellAmpApprox: invalid input or below threshold\n";
#endif
        return {0.0, 0.0};
    }

    // r = 1 - beta^2
    const double s_hat = sqrt_s_hat * sqrt_s_hat;
    const double r = 2.0 * (m1 * m1 + m2 * m2) / s_hat;
    if (r > 1.0) {
#ifdef DEBUG
        std::cerr << "offShellAmpApprox: r > 1.0\n";
#endif
        return {0.0, 0.0};
    }

    const double beta2 = std::max(0.0, 1.0 - r);  // numerical safety

    const double cos_th2 = cos_th * cos_th;
    const double sin_th2 = std::max(0.0, 1.0 - cos_th2);  // numerical safety

    // sin^2 + r*cos^2 is more stable than 1 - beta^2*cos^2
    const double denom = sin_th2 + r * cos_th2;
    if (denom < 1e-18) { return {0.0, 0.0}; }

    const double beta = std::sqrt(beta2);
    const double sqrt_r = std::sqrt(r);  // sqrt(1 - beta2) = sqrt(r)
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
    if (amp == 0.0) {
#ifdef DEBUG
        std::cerr << "offShellAmpApprox: forbidden helicity combinations\n";
#endif
        return {0.0, 0.0};
    }

    return {(COUPLING_FACTOR / denom) * amp, 0.0};
}

PolarizationCoefficients computeCoeffsOnShell(double sqrt_s_hat,
                                              double cos_th) {
    if (sqrt_s_hat < TTBARTHRES) { return {}; }

    return averageHelicities(
        [&](Helicity l1, Helicity l2) -> PolarizationCoefficients {
            const auto pp = onShellHelAmp(sqrt_s_hat, cos_th, l1, l2,
                                          Helicity::PLUS, Helicity::PLUS);
            const double pp2 = std::norm(pp);
            const auto mm = onShellHelAmp(sqrt_s_hat, cos_th, l1, l2,
                                          Helicity::MINUS, Helicity::MINUS);
            const double mm2 = std::norm(mm);
            const auto mp = onShellHelAmp(sqrt_s_hat, cos_th, l1, l2,
                                          Helicity::MINUS, Helicity::PLUS);
            const double mp2 = std::norm(mp);
            const auto pm = onShellHelAmp(sqrt_s_hat, cos_th, l1, l2,
                                          Helicity::PLUS, Helicity::MINUS);
            const double pm2 = std::norm(pm);

            return {
                pp2 + mm2,                                 // c1
                pp2 - mm2,                                 // c2
                mp2 + pm2,                                 // c3
                -mp2 + pm2,                                // c4
                ((pp - mm) * std::conj(mp - pm)).real(),   // c5
                -((pp - mm) * std::conj(mp + pm)).imag(),  // c6
                ((pp + mm) * std::conj(mp + pm)).real(),   // c7
                -((pp + mm) * std::conj(mp - pm)).imag(),  // c8
                ((pp + mm) * std::conj(mp - pm)).real(),   // c9
                -((pp + mm) * std::conj(mp + pm)).imag(),  // c10
                ((pp - mm) * std::conj(mp + pm)).real(),   // c11
                -((pp - mm) * std::conj(mp - pm)).imag(),  // c12
                -2.0 * (pp * std::conj(mm)).real(),        // c13
                2.0 * (pp * std::conj(mm)).imag(),         // c14
                -2.0 * (mp * std::conj(pm)).real(),        // c15
                -2.0 * (mp * std::conj(pm)).imag()         // c16
            };
        });  // return averageHelicities
}
}  // namespace gagatt
