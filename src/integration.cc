#include "integration.h"
#include <gsl/gsl_integration.h>
#include <functional>
#include <memory>

namespace gagatt {
double inner_m2_integrand(double m2, void *params) {
    auto *ps = static_cast<IntegralOffShell *>(params);
    return ps->func(ps->sqrt_s_hat, ps->cos_th, ps->current_m1, m2);
}

double outer_m1_integrand(double m1, void *params) {
    auto *ps = static_cast<IntegralOffShell *>(params);
    ps->current_m1 = m1;

    // m2 <= sqrt(s) - m1
    const double m2_min = ps->m_threshold;
    const double m2_max = ps->sqrt_s_hat - m1;

    if (m2_max < m2_min) { return 0.0; }

    gsl_function F;
    F.function = &inner_m2_integrand;
    F.params = ps;

    double result = 0.0;
    double error = 0.0;
    gsl_integration_qags(&F, m2_min, m2_max, ps->epsabs, ps->epsrel, 1000,
                         ps->inner_ws, &result, &error);

    return result;
}

double integrateOffShellMasses(
    double sqrt_s_hat, double cos_th, double m_min,
    std::function<double(double, double, double, double)> func,
    double epsrel = 1e-4) {
    if (sqrt_s_hat < 2.0 * m_min) { return 0.0; }

    double result = 0.0;
    double error = 0.0;

    return result;
}
}  // namespace gagatt
