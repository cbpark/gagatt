#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include "spin_density.h"

constexpr double SQRTS_MIN = 346.0;  // just above 2 * m_t
constexpr double SQRTS_MAX = 2000.0;

using namespace gagatt;

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "usage: ./bin/sqrt_s_hat_unpol <output.dat> <cos(Theta)> "
                     "N_SQRTS\n"
                  << "  <output.dat>: output file name.\n"
                  << "  <cos(Theta)>: cos(Theta) value\n"
                  << "  N_SQRTS: the size of grid for sqrt(s)\n";
        return EXIT_FAILURE;
    }

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }
    fout << "# (1) cos(theta) (2) sqrt(s_hat) (3) negativity "
            "(4) concurrence (5) horodecki (6) marker\n";

    const double cos_th = std::stod(argv[2]);
    const int N_SQRTS = std::stoi(argv[3]);
    const int N_TOTAL = N_SQRTS;
    const double d_sqrts = (SQRTS_MAX - SQRTS_MIN) / N_SQRTS;
    for (int j = 0; j < N_SQRTS; ++j) {
        const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;

        const auto sdc = SDMatrixCoefficients(sqrt_s_hat, cos_th);
        const auto rho = spinDensityMatrix(sdc);

        fout << std::format(
            "{:.2f}{:>12.4f}{:>12.8f}{:>12.8f}{:>12.8f}{:>10.6f}\n", cos_th,
            sqrt_s_hat, negativity(rho), getConcurrence(rho),
            horodeckiMeasure(sdc), entanglementMarker(sdc));

        if (const int np = j + 1; np % 1000 == 0) {
            std::cout << std::format("-- progress: {} / {}\n", np, N_TOTAL);
        }
    }
    fout << '\n';  // blank line = scanline separator for pm3d

    std::cout << "sqrt_s_hat_unpol: " << argv[1] << " written.\n";
}
