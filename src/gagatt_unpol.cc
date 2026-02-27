#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include "spin_density.h"

constexpr double COS_TH_MIN = -1.0;
constexpr double COS_TH_MAX = 1.0;
constexpr double SQRTS_MIN = 346.0;  // just above 2 * m_t
constexpr double SQRTS_MAX = 2000.0;

using namespace gagatt;

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "usage: ./bin/gagatt_unpol <output.dat> N_COS N_SQRTS\n"
                  << "  <output.dat>: output file name.\n"
                  << "  N_COS: the size of grid for cos(Theta)\n"
                  << "  N_SQRTS: the size of grid for sqrt(s)\n";
        return EXIT_FAILURE;
    }

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }
    fout << "# (1) cos(theta) (2) sqrt(s_hat) (3) rho_t2 (4) concurrence (5) "
            "D (6) bell\n";

    const int N_COS = std::stoi(argv[2]);
    const int N_SQRTS = std::stoi(argv[3]);
    const int N_TOTAL = N_COS * N_SQRTS;
    const double d_cos = (COS_TH_MAX - COS_TH_MIN) / N_COS;
    const double d_sqrts = (SQRTS_MAX - SQRTS_MIN) / N_SQRTS;

    for (int i = 0; i < N_COS; ++i) {
        const double cos_th = COS_TH_MIN + (i + 0.5) * d_cos;

        for (int j = 0; j < N_SQRTS; ++j) {
            const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;

            const auto sdc = SDMatrixCoefficients(sqrt_s_hat, cos_th);
            const auto rho = spinDensityMatrix(sdc);

            fout << std::format(
                "{:.2f}\t{:>10.4f}\t{:d}\t{:d}\t{:d}\t{:d}\n", cos_th, sqrt_s_hat,
                isEntangledByPH(rho), isEntangledByConcurrence(rho),
                isEntangledByD(sdc), violatesBellInequality(sdc));

            if (const int np = i * N_SQRTS + j + 1; np % 1000 == 0) {
                std::cout << std::format("-- progress: {} / {}\n", np, N_TOTAL);
            }
        }
        fout << '\n';  // blank line = scanline separator for pm3d
    }

    std::cout << "gagatt_unpol: the output has been stored in " << argv[1]
              << '\n';
}
