#ifndef SRC_INTEGRATION_H
#define SRC_INTEGRATION_H

#include <gsl/gsl_integration.h>
#include <functional>  // std::function

namespace gagatt {
// Parameters to pass into the integrand
struct IntegralOffShell {
    // The function to integrate: f(sqrt_s, cos_th, m1, m2)
    std::function<double(double, double, double, double)> func;

    double sqrt_s_hat;   // Collision energy
    double cos_th;       // Scattering angle
    double m_threshold;  // Minimum possible mass

    // GSL workspace for the inner integral
    gsl_integration_workspace *inner_ws;

    // Current integration state
    double current_m1 = 0.0;

    // Error control
    double epsabs = 1e-9;
    double epsrel = 1e-6;
};

// double inner_m2_integrand(double m2, void *params);

// double outer_m1_integrand(double m1, void *params);

double integrateOffShellMasses(
    double sqrt_s_hat, double cos_th, double m_min,
    std::function<double(double, double, double, double)> func,
    double epsrel = 1e-6);
}  // namespace gagatt

#endif  // SRC_INTEGRATION_H
