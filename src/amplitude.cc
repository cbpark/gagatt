#include "amplitude.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include "constants.h"
#include "helicity.h"
#include "photon.h"

#ifdef DEBUG
#include <iostream>
#endif

namespace gagatt {
// Internal helper to avoid redundant kinematic calculations
struct KinematicContext {
    double beta, sqrt_r, overall_fac, sin_th, sin_th2, cos_th;
    bool valid;

    KinematicContext(double sqrt_s_hat, double ct, double m1, double m2)
        : cos_th(ct) {
        // phase-space threshold
        if (sqrt_s_hat < m1 + m2 || std::abs(ct) > 1.0) {
            valid = false;
            return;
        }

        const double s_hat = sqrt_s_hat * sqrt_s_hat;
        const double r = 2.0 * (m1 * m1 + m2 * m2) / s_hat;
        // r > 1 can only happen below threshold (already guarded above for the
        // equal-mass case).  Guard for the asymmetric case.
        if (r > 1.0) {
            valid = false;
            return;
        }

        sin_th2 = std::max(0.0, 1.0 - ct * ct);

        const double denom = sin_th2 + r * ct * ct;
        if (denom < 1e-18) {
            valid = false;
            return;
        }

        valid = true;
        beta = std::sqrt(1.0 - r);
        sqrt_r = std::sqrt(r);
        sin_th = std::sqrt(sin_th2);
        overall_fac = COUPLING_FACTOR / denom;
    }
};

// Leading-order helicity amplitude (real at tree-level)
Amplitude computeAmp(const KinematicContext &k, Helicity l1, Helicity l2,
                     Helicity s1, Helicity s2) {
    if (!k.valid) { return {0.0, 0.0}; }

    const double l1_val = toDouble(l1);
    const double s1_val = toDouble(s1);
    double amp = 0.0;

    if (l1 == l2) {
        // identical photon helicities
        if (s1 == s2) {  // identical top helicities
            amp = k.sqrt_r * (k.beta * s1_val + l1_val);
        }
    } else {
        // opposite photon helicities
        if (s1 == s2) {  // identical top helicities
            amp = -k.beta * k.sqrt_r * s1_val * k.sin_th2;
        } else {  // opposite top helicities
            amp = -k.beta * k.sin_th * (s1_val * l1_val + k.cos_th);
        }
    }
    return {k.overall_fac * amp, 0.0};
}

// Per-(l1, l2) polarization coefficients (no averaging, no weighting).
PolarizationCoefficients polCoeffsForHelicity(const KinematicContext &k,
                                              Helicity l1, Helicity l2) {
    if (!k.valid) { return {}; }

    const auto pp = computeAmp(k, l1, l2, Helicity::PLUS, Helicity::PLUS);
    const auto mm = computeAmp(k, l1, l2, Helicity::MINUS, Helicity::MINUS);
    const auto mp = computeAmp(k, l1, l2, Helicity::MINUS, Helicity::PLUS);
    const auto pm = computeAmp(k, l1, l2, Helicity::PLUS, Helicity::MINUS);

    const double pp2 = std::norm(pp), mm2 = std::norm(mm);
    const double mp2 = std::norm(mp), pm2 = std::norm(pm);

    const auto pp_plus_mm = pp + mm;
    const auto pp_minus_mm = pp - mm;
    const auto mp_plus_pm = mp + pm;
    const auto mp_minus_pm = mp - pm;

    return {
        pp2 + mm2,                                       // c1
        pp2 - mm2,                                       // c2
        mp2 + pm2,                                       // c3
        pm2 - mp2,                                       // c4  (= -mp2 + pm2)
        (pp_minus_mm * std::conj(mp_minus_pm)).real(),   // c5
        -(pp_minus_mm * std::conj(mp_plus_pm)).imag(),   // c6
        (pp_plus_mm * std::conj(mp_plus_pm)).real(),     // c7
        -(pp_plus_mm * std::conj(mp_minus_pm)).imag(),   // c8
        (pp_plus_mm * std::conj(mp_minus_pm)).real(),    // c9
        -(pp_plus_mm * std::conj(mp_plus_pm)).imag(),    // c10
        (pp_minus_mm * std::conj(mp_plus_pm)).real(),    // c11
        -(pp_minus_mm * std::conj(mp_minus_pm)).imag(),  // c12
        -2.0 * (pp * std::conj(mm)).real(),              // c13
        2.0 * (pp * std::conj(mm)).imag(),               // c14
        -2.0 * (mp * std::conj(pm)).real(),              // c15
        -2.0 * (mp * std::conj(pm)).imag()               // c16
    };
}

PolarizationCoefficients polCoeffsForHelicity(double sqrt_s_hat, double cos_th,
                                              Helicity l1, Helicity l2) {
    KinematicContext k{sqrt_s_hat, cos_th, MTOP, MTOP};
    return polCoeffsForHelicity(k, l1, l2);
}

Amplitude offShellAmpApprox(double sqrt_s_hat, double cos_th, double m1,
                            double m2, Helicity l1, Helicity l2, Helicity s1,
                            Helicity s2) {
    KinematicContext k(sqrt_s_hat, cos_th, m1, m2);
    return computeAmp(k, l1, l2, s1, s2);
}

// uniform 1/4 average over helicities.
PolarizationCoefficients computePolCoeffs(double sqrt_s_hat, double cos_th) {
    KinematicContext k(sqrt_s_hat, cos_th, MTOP, MTOP);
    return averageHelicities([&](Helicity l1, Helicity l2) {
        return polCoeffsForHelicity(sqrt_s_hat, cos_th, l1, l2);
    });
}

// weighted version.
// index: 0 = (+,+), 1 = (+,−), 2 = (−,+), 3 = (−,−).
PolarizationCoefficients computePolCoeffsWeighted(
    double sqrt_s_hat, double cos_th, const std::array<double, 4> &weights) {
    KinematicContext k(sqrt_s_hat, cos_th, MTOP, MTOP);
    if (!k.valid) { return {}; }

    constexpr std::pair<Helicity, Helicity> hels[4] = {
        {Helicity::PLUS, Helicity::PLUS},
        {Helicity::PLUS, Helicity::MINUS},
        {Helicity::MINUS, Helicity::PLUS},
        {Helicity::MINUS, Helicity::MINUS},
    };

#ifdef DEBUG
    std::cout << "computePolCoeffsWeighted: weight(++, +-, -+, --) = "
              << weights[0] << ", " << weights[1] << ", " << weights[2] << ", "
              << weights[3] << '\n';
#endif

    PolarizationCoefficients total{};
    for (int i = 0; i < 4; ++i) {
        if (weights[i] < 1e-12) { continue; }  // skip negligible weights
        total +=
            polCoeffsForHelicity(k, hels[i].first, hels[i].second) * weights[i];
    }
    return total;
}
}  // namespace gagatt
