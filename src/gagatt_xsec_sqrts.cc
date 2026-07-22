// gagatt_xsec_sqrts.cc
//
// Compute the total cross section sigma(e+e- -> gamma gamma -> t tbar) [fb]
// as a function of sqrt(s), the e+e- CM energy, by scanning over sqrt(s)
// and integrating over sqrt(s_hat) and cos(Theta) at each point.
//
// The total cross section at a given sqrt(s) is:
//   sigma = integral_{2M_t}^{sqrt(s)} d sqrt(s_hat)
//           integral_{-1}^{1} d cos(Theta)
//           sum_{l1,l2} L^{l1,l2}(z) * d sigma_hat^{l1,l2} / d cos(Theta)
//           * GEV2_TO_FB / sqrt(s)
//
// where z = sqrt(s_hat) / sqrt(s).
//
// The effective cross section (for the dileptonic final state) is
// sigma_eff = sigma * BR(t tbar -> l+ l-).
//
// Usage:
//   ./bin/gagatt_xsec_sqrts <output> <pe1> <pe2> <N_sqrts> [N_COS]
//   [N_SQRTS_HAT]

#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>
#include "amplitude.h"
#include "constants.h"
#include "helicity.h"
#include "mc.h"         // X_DEFAULT
#include "mc_helper.h"  // partialXsec, eventRate
#include "photon.h"     // lumiWeightsAndTotal
#include "spin_density.h"

using namespace gagatt;

// ---------------------------------------------------------------------------
// diffEventRate_fixedHel:
//   d^2 sigma^{l1,l2} / (d sqrt(s_hat) d cos(Theta))  [fb/GeV]
// ---------------------------------------------------------------------------
double diffEventRate_fixedHel(double sqrt_s_hat, double cos_th, double L_ll,
                              Helicity l1, Helicity l2, double sqrt_s) {
    if (sqrt_s_hat < TTBARTHRES || L_ll <= 0.0) { return 0.0; }
    const auto pol = polCoeffsForHelicity(sqrt_s_hat, cos_th, l1, l2);
    const SDMatrixCoefficients sdc(pol);
    return eventRate(sqrt_s_hat, sdc, L_ll, sqrt_s) * 4.0;
}

// Integrate over cos_th using Simpson's rule.
template <typename Func>
double integrate_cos(Func f, double cos_min, double cos_max, int n) {
    if (n % 2 == 1) ++n;
    const double d = (cos_max - cos_min) / n;
    double sum = f(cos_min) + f(cos_max);
    for (int i = 1; i < n; ++i) {
        const double c = cos_min + i * d;
        sum += (i % 2 == 0) ? 2.0 * f(c) : 4.0 * f(c);
    }
    return sum * d / 3.0;
}

// ---------------------------------------------------------------------------
// computeTotalXsec:
//   Compute sigma(e+e- -> gamma gamma -> t tbar) for a given sqrt(s).
//
//   Integrates over sqrt(s_hat) in [sqrts_hat_min, sqrt(s)] using the
//   midpoint rule, and over cos(Theta) in [-1, 1] using Simpson's rule.
//   Returns {sigma_total [fb], sigma_++ , sigma_+-, sigma_-+, sigma_--}.
// ---------------------------------------------------------------------------
struct XsecResult {
    double total;
    double pp, pm, mp, mm;
};

XsecResult computeTotalXsec(double sqrt_s, double pe1, double pe2, double x,
                            int n_sqrts_hat, int n_cos) {
    const double pc1 = -pe1;
    const double pc2 = -pe2;
    const double sqrts_hat_min = 2.0 * MTOP + 1.0;
    const double sqrts_hat_max = sqrt_s;
    const double d_sqrts_hat = (sqrts_hat_max - sqrts_hat_min) / n_sqrts_hat;

    XsecResult res = {0.0, 0.0, 0.0, 0.0, 0.0};

    for (int j = 0; j < n_sqrts_hat; ++j) {
        const double sqrt_s_hat = sqrts_hat_min + (j + 0.5) * d_sqrts_hat;
        const double z = sqrt_s_hat / sqrt_s;

        const auto [lw, L_tot] = lumiWeightsAndTotal(z, x, pe1, pc1, pe2, pc2);

        const double L_pp = lw.w[0] * L_tot;
        const double L_pm = lw.w[1] * L_tot;
        const double L_mp = lw.w[2] * L_tot;
        const double L_mm = lw.w[3] * L_tot;

        auto f_pp = [&](double cth) -> double {
            return diffEventRate_fixedHel(sqrt_s_hat, cth, L_pp, Helicity::PLUS,
                                          Helicity::PLUS, sqrt_s);
        };
        auto f_pm = [&](double cth) -> double {
            return diffEventRate_fixedHel(sqrt_s_hat, cth, L_pm, Helicity::PLUS,
                                          Helicity::MINUS, sqrt_s);
        };
        auto f_mp = [&](double cth) -> double {
            return diffEventRate_fixedHel(
                sqrt_s_hat, cth, L_mp, Helicity::MINUS, Helicity::PLUS, sqrt_s);
        };
        auto f_mm = [&](double cth) -> double {
            return diffEventRate_fixedHel(sqrt_s_hat, cth, L_mm,
                                          Helicity::MINUS, Helicity::MINUS,
                                          sqrt_s);
        };

        const double r_pp = integrate_cos(f_pp, -1.0, 1.0, n_cos);
        const double r_pm = integrate_cos(f_pm, -1.0, 1.0, n_cos);
        const double r_mp = integrate_cos(f_mp, -1.0, 1.0, n_cos);
        const double r_mm = integrate_cos(f_mm, -1.0, 1.0, n_cos);

        res.pp += r_pp * d_sqrts_hat;
        res.pm += r_pm * d_sqrts_hat;
        res.mp += r_mp * d_sqrts_hat;
        res.mm += r_mm * d_sqrts_hat;
    }

    res.total = res.pp + res.pm + res.mp + res.mm;
    return res;
}

// sqrt(s) scan range.
// The photon luminosity is nonzero only when z < x/(1+x).
// For x = 4.8, z_max = 4.8/5.8 ~ 0.83, so sqrt(s) must exceed
// 2M_t / z_max ~ 417 GeV for any t tbar production.
constexpr double SQRTS_SCAN_MIN = 350.0;   // GeV
constexpr double SQRTS_SCAN_MAX = 2000.0;  // GeV

int main(int argc, char *argv[]) {
    if (argc != 5 && argc != 6 && argc != 7) {
        std::cerr << "usage: ./bin/gagatt_xsec_sqrts <output> <pe1> <pe2> "
                     "<N_sqrts> [N_COS] [N_SQRTS_HAT]\n"
                  << " <output>      : output file\n"
                  << " <pe1>         : electron polarization\n"
                  << " <pe2>         : positron polarization\n"
                  << " <N_sqrts>     : number of sqrt(s) scan points\n"
                  << " [N_COS]       : number of cos integration steps "
                     "(default 1000)\n"
                  << " [N_SQRTS_HAT] : number of sqrt(s_hat) integration bins "
                     "(default 200)\n";
        return EXIT_FAILURE;
    }

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    const double pe1 = std::stod(argv[2]);
    const double pe2 = std::stod(argv[3]);
    const int N_SQRTS = std::stoi(argv[4]);
    const int N_COS = (argc >= 6) ? std::stoi(argv[5]) : 1000;
    const int N_SQRTS_HAT = (argc >= 7) ? std::stoi(argv[6]) : 200;
    const double x = X_DEFAULT;

    fout << "# sigma(e+e- -> gamma gamma -> t tbar) [fb] vs sqrt(s)\n";
    fout << "# pe1 = " << pe1 << ", pe2 = " << pe2 << ", x = " << x << "\n";
    fout << "# BR(t tbar -> l+ l-) = " << BRLL << "\n";
    fout << "# sqrt(s) in [" << SQRTS_SCAN_MIN << ", " << SQRTS_SCAN_MAX
         << "] GeV\n";
    fout << "# sigma = sigma_++ + sigma_+- + sigma_-+ + sigma_--\n";
    fout << std::format(
        "# {:>12s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s}\n", "sqrt_s",
        "sigma", "sigma_eff", "sigma_++", "sigma_+-", "sigma_-+", "sigma_--");

    const double d_sqrts = (SQRTS_SCAN_MAX - SQRTS_SCAN_MIN) / N_SQRTS;

    for (int k = 0; k < N_SQRTS; ++k) {
        const double sqrt_s = SQRTS_SCAN_MIN + (k + 0.5) * d_sqrts;

        const XsecResult res =
            computeTotalXsec(sqrt_s, pe1, pe2, x, N_SQRTS_HAT, N_COS);

        const double sigma_eff = res.total * BRLL;

        fout << std::format(
            "{:14.2f} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e}\n",
            sqrt_s, res.total, sigma_eff, res.pp, res.pm, res.mp, res.mm);

        if ((k + 1) % 10 == 0 || k == 0) {
            std::cout << std::format(
                "-- {}/{} sqrt_s={:.1f} GeV  sigma={:.6e} fb  "
                "sigma_eff={:.6e} fb\n",
                k + 1, N_SQRTS, sqrt_s, res.total, sigma_eff);
        }
    }

    fout.close();
    std::cout << "gagatt_xsec_sqrts: " << argv[1] << '\n';
    return EXIT_SUCCESS;
}
