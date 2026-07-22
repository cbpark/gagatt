#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include "amplitude.h"
#include "constants.h"
#include "helicity.h"
#include "mc_helper.h"  // partialXsec
#include "spin_density.h"

using namespace gagatt;

double diffXsec_fixedHel(double sqrt_s_hat, double cos_th, Helicity l1,
                         Helicity l2) {
    if (sqrt_s_hat < TTBARTHRES) { return 0.0; }
    // per-(l1,l2) polarisation coefficients
    const auto pol = polCoeffsForHelicity(sqrt_s_hat, cos_th, l1, l2);
    // build SDMatrixCoefficients from pol -> norm = c1+c3
    const SDMatrixCoefficients sdc(pol);
    const double dsig_GeV2 = partialXsec(sqrt_s_hat, sdc);  // GeV^-2
    return dsig_GeV2 * 4.0 * GEV2_TO_FB;                    // fb per unit cos
}

double diffXsec_unpol(double sqrt_s_hat, double cos_th) {
    if (sqrt_s_hat < TTBARTHRES) { return 0.0; }
    // SDMatrixCoefficients(sqrt,cos) internally calls computePolCoeffs()
    // which sums over l1,l2 (not averaged)
    const SDMatrixCoefficients sdc(sqrt_s_hat, cos_th);
    const double dsig_GeV2_summed = partialXsec(sqrt_s_hat, sdc);
    // average over 4 photon helicity combos for unpolarized
    return dsig_GeV2_summed * GEV2_TO_FB;
}

// Integrate over cos_th in [cos_min, cos_max] using Simpson's rule
template <typename Func>
double integrate_cos(Func f, double cos_min, double cos_max, int n) {
    if (n % 2 == 1) ++n;
    const double d = (cos_max - cos_min) / n;
    double sum = f(cos_min) + f(cos_max);
    for (int i = 1; i < n; ++i) {
        const double c = cos_min + i * d;
        const double fi = f(c);
        sum += (i % 2 == 0) ? 2.0 * fi : 4.0 * fi;
    }
    return sum * d / 3.0;
}

constexpr double SQRTS_MIN = 2.0 * MTOP + 1.0;
constexpr double SQRTS_MAX = 1000.0;
constexpr double COS_MIN = -1.0;
constexpr double COS_MAX = 1.0;

int main(int argc, char *argv[]) {
    if (argc != 3 && argc != 4) {
        std::cerr
            << "usage: ./bin/gagatt_xsec_parton <output.dat> <N_SQRTS> "
               "[N_COS]\n"
            << "  <output.dat>: output file\n"
            << "  <N_SQRTS>: number of sqrt(s_hat) points\n"
            << "  [N_COS]: number of cos integration steps (default 1000)\n";
        return EXIT_FAILURE;
    }

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    const int N_SQRTS = std::stoi(argv[2]);
    const int N_COS = (argc == 4) ? std::stoi(argv[3]) : 1000;

    fout << "# sqrt(s_hat) vs partial cross section integrated over "
            "cos(Theta) in ["
         << COS_MIN << ", " << COS_MAX << "]\n";
    fout << "# MTOP = " << MTOP << " GeV, NC = " << NC << ", ALPHA = " << ALPHA
         << "\n";
    fout << "# d sigma / d cos [fb] = partialXsec/4 * GEV2->FB (unpol)\n";
    fout
        << "# sigma_hat = int_{cos_min}^{cos_max} d sigma / d cos d cos [fb]\n";
    fout << std::format("# {:>12s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s}\n",
                        "sqrt_s_hat", "sigma_unpol", "sigma_++", "sigma_+-",
                        "sigma_-+", "sigma_--");

    const double d_sqrts = (SQRTS_MAX - SQRTS_MIN) / N_SQRTS;
    for (int j = 0; j < N_SQRTS; ++j) {
        const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;

        auto f_unpol = [&](double cth) {
            return diffXsec_unpol(sqrt_s_hat, cth);
        };
        auto f_pp = [&](double cth) {
            return diffXsec_fixedHel(sqrt_s_hat, cth, Helicity::PLUS,
                                     Helicity::PLUS);
        };
        auto f_pm = [&](double cth) {
            return diffXsec_fixedHel(sqrt_s_hat, cth, Helicity::PLUS,
                                     Helicity::MINUS);
        };
        auto f_mp = [&](double cth) {
            return diffXsec_fixedHel(sqrt_s_hat, cth, Helicity::MINUS,
                                     Helicity::PLUS);
        };
        auto f_mm = [&](double cth) {
            return diffXsec_fixedHel(sqrt_s_hat, cth, Helicity::MINUS,
                                     Helicity::MINUS);
        };

        const double sig_unpol =
            integrate_cos(f_unpol, COS_MIN, COS_MAX, N_COS);
        const double sig_pp = integrate_cos(f_pp, COS_MIN, COS_MAX, N_COS);
        const double sig_pm = integrate_cos(f_pm, COS_MIN, COS_MAX, N_COS);
        const double sig_mp = integrate_cos(f_mp, COS_MIN, COS_MAX, N_COS);
        const double sig_mm = integrate_cos(f_mm, COS_MIN, COS_MAX, N_COS);

        fout << std::format(
            "{:14.4f} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e}\n",
            sqrt_s_hat, sig_unpol, sig_pp, sig_pm, sig_mp, sig_mm);

        if ((j + 1) % 100 == 0) {
            std::cout << std::format("-- {}/{}  sqrt={:.1f}  sig={:.4e} fb\n",
                                     j + 1, N_SQRTS, sqrt_s_hat, sig_unpol);
        }
    }

    std::cout << "sqrt_s_hat_xsec: " << argv[1] << '\n';
    return EXIT_SUCCESS;
}
