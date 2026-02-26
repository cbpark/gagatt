#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include "spin_density.h"

constexpr int N_COS = 10;
constexpr int N_SQRTS = 10;
constexpr int N_TOTAL = N_COS * N_SQRTS;

constexpr double COS_TH_MIN = -1.0;
constexpr double COS_TH_MAX = 1.0;
constexpr double SQRT_TAU_MIN = 0.0;
constexpr double SQRT_TAU_MAX = 1.0;
constexpr double d_cos = (COS_TH_MAX - COS_TH_MIN) / N_COS;
constexpr double d_sqrt_tau = (SQRT_TAU_MAX - SQRT_TAU_MIN) / N_SQRTS;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "usage: ./bin/gagatt_pol <output.dat> <sqrt_s>\n";
        return EXIT_FAILURE;
    }

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    const double sqrt_s = std::stod(argv[2]);
    std::cout << "gagatt_pol: sqrt(s) = " << sqrt_s << '\n';

    for (int i = 0; i < N_COS; ++i) {
        const double cos_th = COS_TH_MIN + (i + 0.5) * d_cos;

        for (int j = 0; j < N_SQRTS; ++j) {
            // const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;
            const double sqrt_tau = 0.0 + (j + 0.5) * d_sqrt_tau;
            const double sqrt_s_hat = sqrt_s * sqrt_tau;

            fout << std::format("{:.2f}\t{:>12.4f}\n", cos_th, sqrt_s_hat);

            if (const int np = i * N_SQRTS + j + 1; np % 1000 == 0) {
                std::cout << std::format("-- progress: {} / {}\n", np, N_TOTAL);
            }
        }
        fout << '\n';  // blank line = scanline separator for pm3d
    }

    std::cout << "gagatt_pol: " << argv[1] << " written.\n";
}
