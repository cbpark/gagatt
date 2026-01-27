#include "integration.h"

namespace gagatt {
double inner_m2_integrand(double m2, void *params) {
    auto *par = static_cast<IntegralOffShell *>(params);
    return par->func(par->sqrt_s_hat, par->cos_th, par->current_m1, m2);
}

double outer_m1_integrand(double m1, void *params) { return 0.0; }
}  // namespace gagatt
