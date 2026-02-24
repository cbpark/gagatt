#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include "spin_density.h"

constexpr int N_COS = 10;
constexpr int N_SQRTS = 10;
constexpr int N_TOTAL = N_COS * N_SQRTS;

constexpr double COS_TH_MIN = -1.0;
constexpr double COS_TH_MAX = 1.0;
constexpr double SQRTS_MIN = 346.0;  // just above 2 * m_t
constexpr double SQRTS_MAX = 2000.0;
constexpr double d_cos = (COS_TH_MAX - COS_TH_MIN) / N_COS;
constexpr double d_sqrts = (SQRTS_MAX - SQRTS_MIN) / N_SQRTS;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "usage: ./bin/gagatt_pol <output.dat>\n";
        return EXIT_FAILURE;
    }

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    for (int i = 0; i < N_COS; ++i) {
        const double cos_th = COS_TH_MIN + (i + 0.5) * d_cos;

        for (int j = 0; j < N_SQRTS; ++j) {
            const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;

            fout << std::format("{:.2f}\t{:>12.4f}\n", cos_th, sqrt_s_hat);

            if (const int np = i * N_SQRTS + j + 1; np % 1000 == 0) {
                std::cout << std::format("-- progress: {} / {}\n", np, N_TOTAL);
            }
        }
        fout << '\n';  // blank line = scanline separator for pm3d
    }

    std::cout << "gagatt_pol: the output has been stored in " << argv[1] << '\n';
}
