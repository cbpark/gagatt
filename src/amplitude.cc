#include "amplitude.h"
#include <cmath>
#include <complex>
#include <utility>
#include "constants.h"
#include "helicity.h"

namespace gagatt {
// Internal helper to avoid redundant kinematic calculations
struct KinematicContext {
    double beta, sqrt_r, overall_fac, sin_th, sin_th2, cos_th;
    bool valid;

    KinematicContext(double sqrt_s_hat, double ct, double m1, double m2)
        : cos_th(ct) {
        const double threshold = m1 + m2;
        if (sqrt_s_hat < threshold || std::abs(ct) > 1.0) {
            valid = false;
            return;
        }

        const double s_hat = sqrt_s_hat * sqrt_s_hat;
        const double r = 2.0 * (m1 * m1 + m2 * m2) / s_hat;
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

Amplitude computeAmp(const KinematicContext &ctx, Helicity l1, Helicity l2,
                     Helicity s1, Helicity s2) {
    if (!ctx.valid) { return {0.0, 0.0}; }

    const double l1_val = toDouble(l1);
    const double s1_val = toDouble(s1);
    double amp = 0.0;

    if (l1 == l2) {
        if (s1 == s2) { amp = ctx.sqrt_r * (ctx.beta * s1_val + l1_val); }
    } else {
        if (s1 == s2) {
            amp = -ctx.beta * ctx.sqrt_r * s1_val * ctx.sin_th2;
        } else {
            amp = -ctx.beta * ctx.sin_th * (s1_val * l1_val + ctx.cos_th);
        }
    }
    return {ctx.overall_fac * amp, 0.0};
}

// Per-(l1, l2) polarization coefficients — no averaging, no weighting.
PolarizationCoefficients polCoeffsForHelicity(const KinematicContext &ctx,
                                              Helicity l1, Helicity l2) {
    const auto pp = computeAmp(ctx, l1, l2, Helicity::PLUS, Helicity::PLUS);
    const auto mm = computeAmp(ctx, l1, l2, Helicity::MINUS, Helicity::MINUS);
    const auto mp = computeAmp(ctx, l1, l2, Helicity::MINUS, Helicity::PLUS);
    const auto pm = computeAmp(ctx, l1, l2, Helicity::PLUS, Helicity::MINUS);

    const double pp2 = std::norm(pp), mm2 = std::norm(mm);
    const double mp2 = std::norm(mp), pm2 = std::norm(pm);

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
}

// uniform 1/4 average over helicities.
PolarizationCoefficients computePolCoeffs(double sqrt_s_hat, double cos_th) {
    KinematicContext ctx(sqrt_s_hat, cos_th, MTOP, MTOP);
    if (!ctx.valid) { return {}; }

    return averageHelicities([&](Helicity l1, Helicity l2) {
        return polCoeffsForHelicity(ctx, l1, l2);
    });
}

// weighted version.
template <typename W>
PolarizationCoefficients computePolCoeffs(double sqrt_s_hat, double cos_th,
                                          W &&weight) {
    KinematicContext ctx(sqrt_s_hat, cos_th, MTOP, MTOP);
    if (!ctx.valid) { return {}; }

    return weightedHelicities(
        [&](Helicity l1, Helicity l2) {
            return polCoeffsForHelicity(ctx, l1, l2);
        },
        std::forward<W>(weight));
}

// explicit instantiation for the most common case: function pointer.
template PolarizationCoefficients
computePolCoeffs<double (*)(Helicity, Helicity)>(double, double,
                                                 double (*&&)(Helicity,
                                                              Helicity));

Amplitude offShellAmpApprox(double sqrt_s_hat, double cos_th, double m1,
                            double m2, Helicity l1, Helicity l2, Helicity s1,
                            Helicity s2) {
    KinematicContext ctx(sqrt_s_hat, cos_th, m1, m2);
    return computeAmp(ctx, l1, l2, s1, s2);
}
}  // namespace gagatt
