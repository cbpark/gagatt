#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "photon.h"
#include "spin_density.h"

constexpr double X = 4.8;

constexpr int N_COS = 200;
constexpr int N_SQRTS = 4000;
constexpr int N_TOTAL = N_COS * N_SQRTS;

constexpr double COS_TH_MIN = -1.0;
constexpr double COS_TH_MAX = 1.0;
constexpr double SQRTS_MIN = 346.0;
constexpr double SQRTS_MAX = 2000.0;
constexpr double d_cos = (COS_TH_MAX - COS_TH_MIN) / N_COS;
constexpr double d_sqrts = (SQRTS_MAX - SQRTS_MIN) / N_SQRTS;

using namespace gagatt;

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "usage: ./bin/gagatt_pol <output.dat>"
                     " <sqrt_s_ee> <Pe>\n"
                  << "  sqrt_s_ee : e+e- CM energy [GeV]\n"
                  << "  Pe        : |Pe| = |Pebar|, with PePc=PebarPcbar=-1\n";
        return EXIT_FAILURE;
    }

    const double sqrt_s_ee = std::stod(argv[2]);
    const double pe_abs = std::stod(argv[3]);

    // PePc = PebarPcbar = -1
    // => (Pe, Pebar) = (+pe_abs, +pe_abs), Pc = -Pe, Pcbar = -Pebar
    const double pe1 = pe_abs, pc1 = -pe_abs;
    const double pe2 = pe_abs, pc2 = -pe_abs;

    std::cout << std::format(
        "gagatt_pol: sqrt_s_ee={:.1f}, Pe={:+.2f}, x={:.1f}\n", sqrt_s_ee,
        pe_abs, X);

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    fout << "# cos_th  sqrt_s_hat  negativity  concurrence"
            "  bell  delta_3\n";

    std::vector<LumiWeights> lumi_cache(N_SQRTS);
    for (int j = 0; j < N_SQRTS; ++j) {
        const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;
        const double sqrt_tau = sqrt_s_hat / sqrt_s_ee;
        lumi_cache[j] = lumiWeights(sqrt_tau, X, pe1, pc1, pe2, pc2);

        if ((j + 1) % 500 == 0) {
            std::cout << std::format("-- lumi cache: {} / {}\n", j + 1,
                                     N_SQRTS);
        }
    }

    for (int i = 0; i < N_COS; ++i) {
        const double cos_th = COS_TH_MIN + (i + 0.5) * d_cos;

        for (int j = 0; j < N_SQRTS; ++j) {
            const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;
            const auto sdc =
                SDMatrixCoefficients(sqrt_s_hat, cos_th, lumi_cache[j]);
            const auto rho = spinDensityMatrix(sdc);

            fout << std::format(
                "{:.4f}\t{:>10.4f}\t{:>12.8f}\t{:>12.8f}\t{:d}\t{:>12.8f}\n",
                cos_th, sqrt_s_hat, negativity(rho), getConcurrence(rho),
                violatesBellInequality(sdc), sdc.cc.trace() / 3.0);

            if (const int np = i * N_SQRTS + j + 1; np % 10000 == 0) {
                std::cout << std::format("-- progress: {} / {}\n", np, N_TOTAL);
            }
        }
        fout << '\n';
    }
    std::cout << std::format("gagatt_pol: output in {}\n", argv[1]);
}
