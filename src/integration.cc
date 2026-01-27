#include "integration.h"

namespace gagatt {
double inner_m2_integrand(double m2, void *params) {
    auto *p = static_cast<CParams *>(params);
    return p->func(p->sqrt_s_hat, p->cos_th, p->current_m1, m2);
}

double outer_m1_integrand(double m1, void *params) { return 0.0; }
}  // namespace gagatt
