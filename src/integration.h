#ifndef SRC_INTEGRATION_H
#define SRC_INTEGRATION_H

#include <gsl/gsl_integration.h>
#include <functional>

namespace gagatt {
// Parameters to pass into the integrand
struct CParams {
    // The function to integrate: f(sqrt_s, cos_th, m1, m2)
    std::function<double(double, double, double, double)> func;

    double sqrt_s_hat;   // Collision energy
    double cos_th;       // Scattering angle
    double m_threshold;  // Minimum possible mass

    // GSL workspace for the inner integral
    gsl_integration_workspace *inner_ws;

    // Current integration state
    double current_m1;

    // Error control
    double epsabs = 1e-8;
    double epsrel = 1e-5;
};

double inner_m2_integrand(double m2, void *params);

double outer_m1_integrand(double m1, void *params);
}  // namespace gagatt

#endif  // SRC_INTEGRATION_H
