#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include "helicity.h"
#include "photon.h"

constexpr double X = 4.8;
constexpr int NP = 250;

using namespace gagatt;

auto parse_double(const char *str, const char *name = "value")
    -> std::optional<double> {
    try {
        return std::stod(str);
    } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument for " << name << ": " << e.what()
                  << '\n';
        return std::nullopt;
    } catch (const std::out_of_range &e) {
        std::cerr << "Out of range for " << name << ": " << e.what() << '\n';
        return std::nullopt;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "usage: ./bin/luminosity_weight <output.dat> <pe1> <pe2>\n"
                  << "  <output.dat>: output file name.\n"
                  << "  <pe1> <pe2>: pols of initial electrons.\n";
        return EXIT_FAILURE;
    }

    auto pe1_ = parse_double(argv[2], "pe1");
    auto pe2_ = parse_double(argv[3], "pe2");
    if (!pe1_ || !pe2_) { return EXIT_FAILURE; }

    const double pe1 = *pe1_;
    const double pe2 = *pe2_;
    const double pc1 = -pe1;
    const double pc2 = -pe2;

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }
    fout << "# (pe1, pe2) = (" << pe1 << ", " << pe2 << ")\n";
    fout << "# (1) sqrt(tau) (2) w++ (3) w+- (4) w-+ (5) w--\n";

    for (int i = 0; i < NP; ++i) {
        // z = sqrt(tau) = sqrt(s_hat / s)
        const double z = static_cast<double>(i + 1) / NP;
        const auto weights = lumiWeights(z, X, pe1, pc1, pe2, pc2);

        fout << z << "  " << weights.wpp << "  " << weights.wpm << "  "
             << weights.wmp << "  " << weights.wmm << '\n';
    }

    std::cout << "luminosity_weight: output in " << argv[1] << '\n';
}
