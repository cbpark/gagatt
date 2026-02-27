#include "photon.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <cmath>
#include <optional>
#include "helicity.h"

namespace gagatt {
double sigmaC(double x, double pe, double pc) {
    if (x < 1e-12) { return 0.0; }

    double sigma_c0, sigma_c1;

    if (x < 1e-4) {
        // Analytical limits to avoid catastrophic cancellation
        sigma_c0 = x * (4.0 / 3.0 - x * (4.0 / 3.0 - 26.0 / 15.0 * x));
        sigma_c1 = x * x * (-1.0 / 3.0 + 5.0 / 6.0 * x);
    } else {
        const double inv_x = 1.0 / x;
        const double inv_x2 = inv_x * inv_x;
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

double fLumi(double x, double y, double pe, double pc) {
    const double y_max = x / (1.0 + x);
    if (y < 0.0 || y >= y_max || y >= 1.0) { return 0.0; }

    const double s_c = sigmaC(x, pe, pc);
    if (s_c <= 0.0) { return 0.0; }

    const double om_y = 1.0 - y;
    const double r = y / (x * om_y);

    const double f = (1.0 / om_y) + om_y - 4.0 * r * (1.0 - r) -
                     r * x * (2.0 * r - 1.0) * (2.0 - y) * pe * pc;

    return f / s_c;
}

// Internal structure to hold pre-calculated constants for integration
struct IntegrationParams {
    double tau;
    double x;
    double pe1, pc1;
    double pe2, pc2;
};

// Internal helper to get both c0 and c2 components for a single photon
struct PhotonComponents {
    double c0, c2;
};

PhotonComponents computePhotonComponents(double x, double y, double pe,
                                         double pc) {
    const double om_y = 1.0 - y;
    const double r = y / (x * om_y);
    const double two_r_1 = 2.0 * r - 1.0;

    double c0 = (1.0 / om_y) + om_y - 4.0 * r * (1.0 - r) -
                r * x * two_r_1 * (2.0 - y) * (pe * pc);

    double c2 = r * x * (1.0 + om_y * two_r_1 * two_r_1) * pe -
                two_r_1 * (1.0 / om_y + om_y) * pc;

    return {c0, c2};
}

// Specific Integrands for Polarization Correlation Terms
double c00Func(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);
    auto p1 = computePhotonComponents(p->x, y, p->pe1, p->pc1);
    auto p2 = computePhotonComponents(p->x, p->tau / y, p->pe2, p->pc2);
    return (p1.c0 * p2.c0) / y;
}

double c20Func(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);
    auto p1 = computePhotonComponents(p->x, y, p->pe1, p->pc1);
    auto p2 = computePhotonComponents(p->x, p->tau / y, p->pe2, p->pc2);
    return (p1.c2 * p2.c0) / y;
}

double c02Func(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);
    auto p1 = computePhotonComponents(p->x, y, p->pe1, p->pc1);
    auto p2 = computePhotonComponents(p->x, p->tau / y, p->pe2, p->pc2);
    return (p1.c0 * p2.c2) / y;
}

double c22Func(double y, void *params) {
    auto *p = static_cast<IntegrationParams *>(params);
    auto p1 = computePhotonComponents(p->x, y, p->pe1, p->pc1);
    auto p2 = computePhotonComponents(p->x, p->tau / y, p->pe2, p->pc2);
    return (p1.c2 * p2.c2) / y;
}

struct PolarizationCorrelation {
    double c00, c20, c02, c22;
    PolarizationCorrelation operator*(double f) const {
        return {c00 * f, c20 * f, c02 * f, c22 * f};
    }
};

PolarizationCorrelation computePolCor(double z, double x, double pe1,
                                      double pc1, double pe2, double pc2) {
    const double ymax = x / (1.0 + x);
    const double tau = z * z;
    const double ymin = tau / ymax;
    if (z <= 0.0 || z >= ymax || ymin >= ymax) { return {}; }

    IntegrationParams params{
        .tau = tau, .x = x, .pe1 = pe1, .pc1 = pc1, .pe2 = pe2, .pc2 = pc2};
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.params = &params;

    auto integrate = [&](auto func_ptr) {
        double res = 0.0, err = 0.0;
        F.function = func_ptr;
        gsl_integration_qags(&F, ymin, ymax, 0, 1e-7, 1000, w, &res, &err);
        return res;
    };

    double c00 = integrate(&c00Func);
    double c20 = integrate(&c20Func);
    double c02 = integrate(&c02Func);
    double c22 = integrate(&c22Func);

    gsl_integration_workspace_free(w);

    return PolarizationCorrelation{c00, c20, c02, c22} * (2.0 * z);
}

double photonLuminosity(double z, double x, double pe1, double pc1, double pe2,
                        double pc2, std::optional<Helicity> l1,
                        std::optional<Helicity> l2) {
    const double sigma_c1 = sigmaC(x, pe1, pc1);
    const double sigma_c2 = sigmaC(x, pe2, pc2);

    if (sigma_c1 <= 0.0 || sigma_c2 <= 0.0) { return 0.0; }

    const auto pol_corr = computePolCor(z, x, pe1, pc1, pe2, pc2);

    const double l1_val = l1 ? toDouble(*l1) : 0.0;
    const double l2_val = l2 ? toDouble(*l2) : 0.0;

    // Combine components: L = (C00 + l1*C20 + l2*C02 + l1*l2*C22) / (sigma1 *
    // sigma2)
    return (pol_corr.c00 + l1_val * pol_corr.c20 + l2_val * pol_corr.c02 +
            l1_val * l2_val * pol_corr.c22) /
           (sigma_c1 * sigma_c2);
}

LumiWeights lumiWeights(double z, double x, double pe1, double pc1, double pe2,
                        double pc2) {
    const double sigma_c1 = sigmaC(x, pe1, pc1);
    const double sigma_c2 = sigmaC(x, pe2, pc2);
    if (sigma_c1 <= 0.0 || sigma_c2 <= 0.0) { return {}; }

    const auto pc = computePolCor(z, x, pe1, pc1, pe2, pc2);
    const double inv_sigma = 1.0 / (sigma_c1 * sigma_c2);

    // L_{l1,l2} = (c00 + l1*c20 + l2*c02 + l1*l2*c22) / (sigma1*sigma2)
    // Order: ++, +-, -+, --
    std::array<double, 4> lumi{};
    int idx = 0;
    for (double l1 : {+1.0, -1.0}) {
        for (double l2 : {+1.0, -1.0}) {
            lumi[idx++] =
                (pc.c00 + l1 * pc.c20 + l2 * pc.c02 + l1 * l2 * pc.c22) *
                inv_sigma;
        }
    }

    // Normalise: w_i = L_i / sum(L)
    const double sum = lumi[0] + lumi[1] + lumi[2] + lumi[3];
    if (sum < 1e-15) { return {}; }

    const double inv_sum = 1.0 / sum;
    return LumiWeights{{lumi[0] * inv_sum, lumi[1] * inv_sum, lumi[2] * inv_sum,
                        lumi[3] * inv_sum}};
}
}  // namespace gagatt
