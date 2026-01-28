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
    std::function<double(double, double, double, double)> func, double epsrel) {
    if (sqrt_s_hat < 2.0 * m_min) { return 0.0; }

    auto outer_ws = std::unique_ptr<gsl_integration_workspace,
                                    void (*)(gsl_integration_workspace *)>(
        gsl_integration_workspace_alloc(2000), gsl_integration_workspace_free);

    auto inner_ws = std::unique_ptr<gsl_integration_workspace,
                                    void (*)(gsl_integration_workspace *)>(
        gsl_integration_workspace_alloc(2000), gsl_integration_workspace_free);

    IntegralOffShell params{.func = std::move(func),
                            .sqrt_s_hat = sqrt_s_hat,
                            .cos_th = cos_th,
                            .m_threshold = m_min,
                            .inner_ws = inner_ws.get(),
                            .epsrel = epsrel};

    gsl_function F;
    F.function = &outer_m1_integrand;
    F.params = &params;

    const double m1_min = m_min;
    const double m1_max = sqrt_s_hat - m_min;

    double result = 0.0;
    double error = 0.0;
    gsl_integration_qags(&F, m1_min, m1_max, 0, epsrel, 2000, outer_ws.get(),
                         &result, &error);

    return result;
}
}  // namespace gagatt
