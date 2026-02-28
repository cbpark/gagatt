#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "photon.h"
#include "spin_density.h"

constexpr double X = 4.8;

constexpr double COS_TH_MIN = -1.0;
constexpr double COS_TH_MAX = 1.0;
// constexpr double SQRTS_MIN = 400.0;
// constexpr double SQRTS_MAX = 400.0;
constexpr double SQRTS_MIN = 346.0;
constexpr double SQRTS_MAX = 2000.0;

using namespace gagatt;

int main(int argc, char *argv[]) {
    if (argc != 7) {
        std::cerr << "usage: ./bin/gagatt_pol <output.dat> N_COS N_SQRTS "
                     "<sqrt_s> <Pe1> <Pe2>\n"
                  << "  N_COS: the size of grid for cos(Theta)\n"
                  << "  N_SQRTS: the size of grid for sqrt(s)\n"
                  << "  sqrt_s : e+e- CM energy [GeV]\n"
                  << "  Pe1, Pe2: helicities of electrons\n";
        return EXIT_FAILURE;
    }

    const double sqrt_s = std::stod(argv[4]);
    const double pe1 = std::stod(argv[5]);
    const double pe2 = std::stod(argv[6]);

    // PePc = PebarPcbar = -1
    const double pc1 = -pe1;
    const double pc2 = -pe2;

    std::cout << std::format(
        "gagatt_pol: sqrt_s={:.1f}, pe1={:.2f}, pe2={:.2f}\n", sqrt_s, pe1,
        pe2);

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }
    fout << "# (1) cos_th  (2) sqrt_s_hat  (3) negativity (4) concurrence (5) "
            "horodecki (6) marker\n";

    const int N_COS = std::stoi(argv[2]);
    const int N_SQRTS = std::stoi(argv[3]);
    const double d_cos = (COS_TH_MAX - COS_TH_MIN) / N_COS;
    const double d_sqrts = (SQRTS_MAX - SQRTS_MIN) / N_SQRTS;
    const int N_TOTAL = N_COS * N_SQRTS;

    std::vector<LumiWeights> lumi_cache(N_SQRTS);
    for (int j = 0; j < N_SQRTS; ++j) {
        const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;
        // std::cout << "sqrt_s_hat: " << sqrt_s_hat << '\n';
        if (sqrt_s_hat > sqrt_s) { break; }

        const double sqrt_tau = sqrt_s_hat / sqrt_s;
        lumi_cache[j] = lumiWeights(sqrt_tau, X, pe1, pc1, pe2, pc2);
        // std::cout << "lumi_cache: " << lumi_cache[j].w[0] << ", "
        //           << lumi_cache[j].w[1] << ", " << lumi_cache[j].w[2] << ", "
        //           << lumi_cache[j].w[3] << '\n';

        if ((j + 1) % 500 == 0) {
            std::cout << std::format("-- lumi cache: {} / {}\n", j + 1,
                                     N_SQRTS);
        }
    }

    for (int i = 0; i < N_COS; ++i) {
        const double cos_th = COS_TH_MIN + (i + 0.5) * d_cos;

        for (int j = 0; j < N_SQRTS; ++j) {
            const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;
            if (sqrt_s_hat > sqrt_s) { break; }

            const auto sdc =
                SDMatrixCoefficients(sqrt_s_hat, cos_th, lumi_cache[j]);
            const auto rho = spinDensityMatrix(sdc);

            fout << std::format(
                "{:.2f}{:>12.4f}{:>12.8f}{:>12.8f}{:>12.8f}{:>10.6f}\n", cos_th,
                sqrt_s_hat, negativity(rho), getConcurrence(rho),
                horodeckiMeasure(sdc), entanglementMarker(sdc));

            if (const int np = i * N_SQRTS + j + 1; np % 1000 == 0) {
                std::cout << std::format("-- progress: {} / {}\n", np, N_TOTAL);
            }
        }
        fout << '\n';
    }
    std::cout << std::format("gagatt_pol: output in {}\n", argv[1]);
}
