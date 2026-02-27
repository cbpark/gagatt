#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include "photon.h"
#include "spin_density.h"

constexpr double COS_TH_MIN = -1.0;
constexpr double COS_TH_MAX = 1.0;
constexpr double SQRTS_MIN = 346.0;  // just above 2 * m_t
constexpr double SQRTS_MAX = 2000.0;

using namespace gagatt;

struct HelicityMode {
    const char *tag;  // "pp", "pm", "mp", "mm"
    Helicity h1, h2;
};

constexpr HelicityMode MODES[] = {{"pp", Helicity::PLUS, Helicity::PLUS},
                                  {"pm", Helicity::PLUS, Helicity::MINUS},
                                  {"mp", Helicity::MINUS, Helicity::PLUS},
                                  {"mm", Helicity::MINUS, Helicity::MINUS}};

LumiWeights weightFor(Helicity h1, Helicity h2) {
    using H = Helicity;
    if (h1 == H::PLUS && h2 == H::PLUS) {
        return {1.0, 0.0, 0.0, 0.0};
    } else if (h1 == H::PLUS && h2 == H::MINUS) {
        return {0.0, 1.0, 0.0, 0.0};
    } else if (h1 == H::MINUS && h2 == H::PLUS) {
        return {0.0, 0.0, 1.0, 0.0};
    } else if (h1 == H::MINUS && h2 == H::MINUS) {
        return {0.0, 0.0, 0.0, 1.0};
    } else {
        return {0.0, 0.0, 0.0, 0.0};
    }
}

void run(const HelicityMode &mode, const std::string &fname, int N_COS,
         int N_SQRTS) {
    std::ofstream fout(fname);
    if (!fout) {
        std::cerr << "Failed to open " << fname << '\n';
        return;
    }

    fout << std::format(
        "# helicity: ({}{})\n"
        "# (1) cos(theta) (2) sqrt(s_hat) (3) negativity "
        "(4) concurrence (5) horodecki (6) marker\n",
        (mode.h1 == Helicity::PLUS ? '+' : '-'),
        (mode.h2 == Helicity::PLUS ? '+' : '-'));

    const int N_TOTAL = N_COS * N_SQRTS;
    const double d_cos = (COS_TH_MAX - COS_TH_MIN) / N_COS;
    const double d_sqrts = (SQRTS_MAX - SQRTS_MIN) / N_SQRTS;

    auto lw = weightFor(mode.h1, mode.h2);

    for (int i = 0; i < N_COS; ++i) {
        const double cos_th = COS_TH_MIN + (i + 0.5) * d_cos;

        for (int j = 0; j < N_SQRTS; ++j) {
            const double sqrt_s_hat = SQRTS_MIN + (j + 0.5) * d_sqrts;
            SDMatrixCoefficients sdc(sqrt_s_hat, cos_th, lw);
            auto rho = spinDensityMatrix(sdc);

            fout << std::format(
                "{:.2f}{:>12.4f}{:>12.8f}{:>12.8f}{:>12.8f}{:>10.6f}\n", cos_th,
                sqrt_s_hat, negativity(rho), getConcurrence(rho),
                horodeckiMeasure(sdc), entanglementMarker(sdc));

            if (const int np = i * N_SQRTS + j + 1; np % 1000 == 0) {
                std::cout << std::format("-- [{}] progress: {} / {}\n",
                                         mode.tag, np, N_TOTAL);
            }
        }
        fout << '\n';  // blank line = scanline separator for pm3d
    }

    std::cout << "gagatt_pol_perfect: " << fname << " written.\n";
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        std::cerr << "usage: ./bin/gagatt_pol_perfect <output.dat> N_COS "
                     "N_SQRTS <mode>\n"
                  << "  <output.dat>: output file name.\n"
                  << "  N_COS: the size of grid for cos(Theta)\n"
                  << "  N_SQRTS: the size of grid for sqrt(s)\n"
                  << "  <mode>: pp | pm | mp | mm\n";
        return EXIT_FAILURE;
    }

    const int N_COS = std::stoi(argv[2]);
    const int N_SQRTS = std::stoi(argv[3]);
    const std::string mode = argv[4];

    bool found = false;
    for (const auto &m : MODES) {
        if (mode == m.tag) {
            run(m, argv[1], N_COS, N_SQRTS);
            found = true;
            break;
        }
    }
    if (!found) {
        std::cerr << "Unknown mode '" << mode
                  << "'. Choose from: pp pm mp mm\n";
        return EXIT_FAILURE;
    }
}
