#include "photon.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <cmath>
#include "amplitude.h"

namespace gagatt {
// (2 pi alpha^2 / (x m_e^2))^{-1} * sigma_c
double sigmaC(double x, double pe, double pc) {
    if (x < 1e-12) return 0.0;

    double sigma_c0, sigma_c1;

    // Taylor expansion for small x to avoid 1/x^2 catastrophic cancellation
    // Analytical limits as x -> 0:
    // sigma_c0 -> (4/3)x - (4/3)x^2 + (26/15)x^3 + ...
    // sigma_c1 -> -(1/3)x^2 + (5/6)x^3 ...
    if (x < 1e-4) {
        sigma_c0 = x * (4.0 / 3.0 - x * (4.0 / 3.0 - 26.0 / 15.0 * x));
        sigma_c1 = x * x * (-1.0 / 3.0 + 5.0 / 6.0 * x);
    } else {
        const double x2 = x * x;
        const double inv_x = 1.0 / x;
        const double inv_x2 = 1.0 / x2;
        const double log1px = std::log1p(x);
        const double inv_1px = 1.0 / (1.0 + x);
        const double inv_1px2 = inv_1px * inv_1px;

        sigma_c0 = (1.0 - 4.0 * inv_x - 8.0 * inv_x2) * log1px +
                   (0.5 + 8.0 * inv_x - 0.5 * inv_1px2);

        sigma_c1 =
            (1.0 + 2.0 * inv_x) * log1px - (2.5 - inv_1px + 0.5 * inv_1px2);
    }

    return sigma_c0 + pe * pc * sigma_c1;
}

// (2 pi alpha^2 / (sigma_c x m_e^2))^{-1} f(y)
double fLumi(double x, double y, double pe, double pc) {
    const double y_max = x / (1.0 + x);
    if (y < 0.0 || y >= y_max || y >= 1.0) { return 0.0; }

    const double sigma_c = sigmaC(x, pe, pc);
    if (sigma_c <= 0.0) { return 0.0; }

    const double om_y = 1.0 - y;
    const double r = y / (x * om_y);

    double f = (1.0 / om_y) + om_y - 4.0 * r * (1.0 - r) -
               r * x * (2.0 * r - 1.0) * (2.0 - y) * pe * pc;

    return f / sigma_c;
}

struct IntegrationParams {
    double tau;
    double x;
    double pe1, pc1;
    double pe2, pc2;
};

double photonLuminosityUnpolFunc(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);

    const double f1 = fLumi(p->x, y, p->pe1, p->pc1);
    if (f1 <= 0.0) { return 0.0; }

    const double f2 = fLumi(p->x, p->tau / y, p->pe2, p->pc2);
    if (f2 <= 0.0) { return 0.0; }

    return (f1 * f2) / y;
}

double photonLuminosityUnpol(double z, double x, double pe1, double pc1,
                             double pe2, double pc2) {
    const double ymax = x / (1.0 + x);
    const double tau = z * z;
    const double ymin = tau / ymax;

    if (z <= 0.0 || z >= ymax || ymin >= ymax) { return 0.0; }

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    IntegrationParams params{
        .tau = tau, .x = x, .pe1 = pe1, .pc1 = pc1, .pe2 = pe2, .pc2 = pc2};

    gsl_function F;
    F.function = &photonLuminosityUnpolFunc;
    F.params = &params;

    double result = 0.0, error = 0.0;
    int status =
        gsl_integration_qags(&F, ymin, ymax, 0, 1e-7, 1000, w, &result, &error);

    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) { return 0.0; }

    // Apply the Jacobian (d_tau / d_z = 2z) here once.
    // Result of integral is dL/d_tau. We return dL/d_z.
    return 2.0 * z * result;
}

double c0(double x, double y, double pe, double pc) {
    const double sigma_c = sigmaC(x, pe, pc);
    return fLumi(x, y, pe, pc) * sigma_c;
}

double c2(double x, double y, double pe, double pc) {
    const double om_y = 1.0 - y;
    const double r = y / (x * om_y);
    const double two_r_1 = 2.0 * r - 1.0;
    return r * x * (1.0 + om_y * two_r_1 * two_r_1) * pe -
           two_r_1 * (1.0 / om_y + om_y) * pc;
}

double c00Func(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);

    const double c_1 = c0(p->x, y, p->pe1, p->pc1);
    if (c_1 <= 0.0) { return 0.0; }

    const double c_2 = c0(p->x, p->tau / y, p->pe2, p->pc2);
    if (c_2 <= 0.0) { return 0.0; }

    return (c_1 * c_2) / y;
}

double c20Func(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);

    const double c_1 = c2(p->x, y, p->pe1, p->pc1);
    if (c_1 <= 0.0) { return 0.0; }

    const double c_2 = c0(p->x, p->tau / y, p->pe2, p->pc2);
    if (c_2 <= 0.0) { return 0.0; }

    return (c_1 * c_2) / y;
}

double c02Func(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);

    const double c_1 = c0(p->x, y, p->pe1, p->pc1);
    if (c_1 <= 0.0) { return 0.0; }

    const double c_2 = c2(p->x, p->tau / y, p->pe2, p->pc2);
    if (c_2 <= 0.0) { return 0.0; }

    return (c_1 * c_2) / y;
}

double c22Func(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);

    const double c_1 = c2(p->x, y, p->pe1, p->pc1);
    if (c_1 <= 0.0) { return 0.0; }

    const double c_2 = c2(p->x, p->tau / y, p->pe2, p->pc2);
    if (c_2 <= 0.0) { return 0.0; }

    return (c_1 * c_2) / y;
}

struct PolarizationCorrelation {
    double c00, c20, c02, c22;

    PolarizationCorrelation operator*(double factor) const {
        return {c00 * factor, c20 * factor, c02 * factor, c22 * factor};
    }
};

PolarizationCorrelation computePolCor(double z, double x, double pe1,
                                      double pc1, double pe2, double pc2) {
    const double ymax = x / (1.0 + x);
    const double tau = z * z;
    const double ymin = tau / ymax;

    if (z <= 0.0 || z >= ymax || ymin >= ymax) { return {}; }

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    IntegrationParams params{
        .tau = tau, .x = x, .pe1 = pe1, .pc1 = pc1, .pe2 = pe2, .pc2 = pc2};

    gsl_function F;
    F.params = &params;

    F.function = &c00Func;
    double c00_res = 0.0, c00_err = 0.0;
    gsl_integration_qags(&F, ymin, ymax, 0, 1e-7, 1000, w, &c00_res, &c00_err);

    F.function = &c20Func;
    double c20_res = 0.0, c20_err = 0.0;
    gsl_integration_qags(&F, ymin, ymax, 0, 1e-7, 1000, w, &c20_res, &c20_err);

    F.function = &c02Func;
    double c02_res = 0.0, c02_err = 0.0;
    gsl_integration_qags(&F, ymin, ymax, 0, 1e-7, 1000, w, &c02_res, &c02_err);

    F.function = &c22Func;
    double c22_res = 0.0, c22_err = 0.0;
    gsl_integration_qags(&F, ymin, ymax, 0, 1e-7, 1000, w, &c22_res, &c22_err);

    gsl_integration_workspace_free(w);

    return {2.0 * z * c00_res, 2.0 * z * c20_res, 2.0 * z * c02_res,
            2.0 * z * c22_res};
}

double photonLuminosity(double z, double x, Helicity l1, Helicity l2) {
    const double pe1 = 1.0, pc1 = -1.0;
    const double pe2 = -1.0, pc2 = 1.0;

    const double sigma_c1 = sigmaC(x, pe1, pc1);
    const double sigma_c2 = sigmaC(x, pe2, pc2);

    const auto pol_corr = computePolCor(z, x, pe1, pc1, pe2, pc2);

    const double l1_val = toDouble(l1);
    const double l2_val = toDouble(l2);

    return (pol_corr.c00 + l1_val * pol_corr.c20 + l2_val * pol_corr.c02 +
            l1_val * l2_val * pol_corr.c22) /
           (sigma_c1 * sigma_c2);
}
}  // namespace gagatt
