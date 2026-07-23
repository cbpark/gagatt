// gagatt_eventrate.cc
//
// Generate the distribution of d sigma / d sqrt(s_hat) [fb/GeV] as a function
// of sqrt(s_hat), integrating over cos(Theta) in [-1, 1], with photon
// luminosity folding.
//
// The differential event rate is:
//   d^2 sigma / (d sqrt(s_hat) d cos(Theta))
//     = sum_{l1,l2} L^{l1,l2}(z) * (d sigma_hat^{l1,l2} / d cos(Theta))
//       * GEV2_TO_FB / sqrt(s)
//
// where z = sqrt(s_hat) / sqrt(s), and L^{l1,l2} is the polarized photon
// luminosity per unit z.
//
// The total event rate (column 2) equals the sum of individual helicity
// contributions (columns 3-6), i.e.,
//   dsig/dsqrts = dsig/dsqrts_++ + dsig/dsqrts_+-
//               + dsig/dsqrts_-+ + dsig/dsqrts_--

#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>
#include "amplitude.h"
#include "constants.h"
#include "helicity.h"
#include "mc_helper.h"  // partialXsec, eventRate
#include "photon.h"     // lumiWeightsAndTotal
#include "spin_density.h"

using namespace gagatt;

constexpr double SQRTS_MIN = 2.0 * MTOP + 1.0;
constexpr double COS_MIN = -1.0;
constexpr double COS_MAX = 1.0;

// ---------------------------------------------------------------------------
// diffEventRate_fixedHel:
//   d^2 sigma^{l1,l2} / (d sqrt(s_hat) d cos(Theta))  [fb/GeV]
//
//   = partialXsec_{l1,l2} * L^{l1,l2} * GEV2_TO_FB / sqrt(s)
//
// This is the contribution of a single photon helicity pair (l1, l2) to the
// differential event rate.  No helicity averaging is applied.
// ---------------------------------------------------------------------------
double diffEventRate_fixedHel(double sqrt_s_hat, double cos_th, double L_ll,
                              Helicity l1, Helicity l2, double sqrt_s) {
    if (sqrt_s_hat < TTBARTHRES || L_ll <= 0.0) { return 0.0; }
    const auto pol = polCoeffsForHelicity(sqrt_s_hat, cos_th, l1, l2);
    const SDMatrixCoefficients sdc(pol);
    return eventRate(sqrt_s_hat, sdc, L_ll, sqrt_s) * 4.0;
}

int main(int argc, char *argv[]) {
    if (argc != 6 && argc != 7) {
        std::cerr
            << "usage: ./bin/gagatt_eventrate <output> <sqrt_s> <pe1> <pe2> "
               "<N_sqrts> [N_COS]\n"
            << " <output> : output file\n"
            << " <sqrt_s> : e+e- CM energy [GeV]\n"
            << " <pe1>    : electron polarization\n"
            << " <pe2>    : positron polarization\n"
            << " <N_sqrts>: number of sqrt(s_hat) points\n"
            << " [N_COS]  : number of cos integration steps (default 1000)\n";
        return EXIT_FAILURE;
    }

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    const double sqrt_s = std::stod(argv[2]);
    const double pe1 = std::stod(argv[3]);
    const double pe2 = std::stod(argv[4]);
    const int N_SQRTS = std::stoi(argv[5]);
    const int N_COS = (argc == 7) ? std::stoi(argv[6]) : 1000;

    const double pc1 = -pe1;
    const double pc2 = -pe2;
    const double sqrts_max = sqrt_s;
    const double x = X_DEFAULT;

    fout << "# d sigma / d sqrt(s_hat) [fb/GeV] vs sqrt(s_hat)\n";
    fout << "# with photon luminosity folding\n";
    fout << "# sqrt_s = " << sqrt_s << " GeV, pe1 = " << pe1
         << ", pe2 = " << pe2 << ", x = " << x << "\n";
    fout << "# sqrt(s_hat) in [" << SQRTS_MIN << ", " << sqrts_max
         << "], cos(Theta) in [" << COS_MIN << ", " << COS_MAX << "]\n";
    fout << "# dsig/dsqrts = sum_{l1,l2} dsig/dsqrts_{l1,l2}\n";
    fout << std::format("# {:>12s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s}\n",
                        "sqrt_s_hat", "dsig/dsqrts", "dsig/dsqrts_++",
                        "dsig/dsqrts_+-", "dsig/dsqrts_-+", "dsig/dsqrts_--");

    const double d_sqrts = (sqrts_max - SQRTS_MIN) / N_SQRTS;
    double total_xsec = 0.0;  // integrated cross section for verification

    for (int j = 0; j < N_SQRTS; ++j) {
        const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;
        const double z = sqrt_s_hat / sqrt_s;

        // Compute luminosity weights and total
        const auto [lw, L_tot] = lumiWeightsAndTotal(z, x, pe1, pc1, pe2, pc2);

        // Individual helicity luminosities: L^{l1,l2} = w[i] * L_tot
        const double L_pp = lw.w[0] * L_tot;
        const double L_pm = lw.w[1] * L_tot;
        const double L_mp = lw.w[2] * L_tot;
        const double L_mm = lw.w[3] * L_tot;

        // Total differential event rate = sum of 4 helicity contributions.
        // We compute each contribution separately and add them, so that
        // dsig/dsqrts = dsig/dsqrts_++ + dsig/dsqrts_+- + dsig/dsqrts_-+
        //             + dsig/dsqrts_--
        // exactly (no factor-of-4 ambiguity from the luminosity-weighted
        // SDMatrixCoefficients normalisation).
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

        const double rate_pp = integrate_cos(f_pp, COS_MIN, COS_MAX, N_COS);
        const double rate_pm = integrate_cos(f_pm, COS_MIN, COS_MAX, N_COS);
        const double rate_mp = integrate_cos(f_mp, COS_MIN, COS_MAX, N_COS);
        const double rate_mm = integrate_cos(f_mm, COS_MIN, COS_MAX, N_COS);

        const double rate_total = rate_pp + rate_pm + rate_mp + rate_mm;
        total_xsec += rate_total * d_sqrts;

        fout << std::format(
            "{:14.4f} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e}\n",
            sqrt_s_hat, rate_total, rate_pp, rate_pm, rate_mp, rate_mm);

        if ((j + 1) % 100 == 0) {
            std::cout << std::format(
                "-- {}/{} sqrts={:.1f} rate={:.4e} fb/GeV\n", j + 1, N_SQRTS,
                sqrt_s_hat, rate_total);
        }
    }

    fout.close();
    std::cout << std::format(
        "total cross section (midpoint rule in sqrt_s_hat) = {:.6f} fb\n",
        total_xsec);
    std::cout << "gagatt_eventrate: " << argv[1] << '\n';
    return EXIT_SUCCESS;
}
