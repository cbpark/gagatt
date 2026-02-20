#include <stdlib.h>
#include <array>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include "helicity.h"
#include "photon.h"

const double X = 4.8;
const int NP = 250;

auto parse_double(const char *str, const char *name = "value")
    -> std::optional<double> {
    try {
        return std::stod(str);
    } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument for " << name << ": " << e.what()
                  << '\n';
        return std::nullopt;
    }
}

void printInput(double pe1, double pe2, double l1, double l2) {
    std::cout << "luminosity: (pe1, pe2) = (" << pe1 << ", " << pe2 << "), "
              << "(l1, l2) = (" << l1 << ", " << l2 << ")\n";
}

int main(int argc, char *argv[]) {
    if (!(argc == 6)) {
        std::cerr
            << "usage: ./bin/luminosity <output.dat> <pe1> <pe2> <l1> <l2>\n"
            << "  <output.dat>: output file name.\n"
            << "  <pe1> <pe2>: pols of initial electrons.\n"
            << "  <l1> <l2>: pols of colliding photons.\n";
        return 1;
    }
    std::ofstream fout(argv[1]);

    auto pe1_ = parse_double(argv[2], "pe1");
    auto pe2_ = parse_double(argv[3], "pe2");
    auto l1_ = parse_double(argv[4], "l1");
    auto l2_ = parse_double(argv[5], "l2");
    if (!pe1_ || !pe2_ || !l1_ || !l2_) { return EXIT_FAILURE; }

    double pe1 = *pe1_;
    double pe2 = *pe2_;
    double pc1 = -pe1;
    double pc2 = -pe2;

    auto l1 = gagatt::toHelicity(l1_.value());
    auto l2 = gagatt::toHelicity(l2_.value());

    printInput(pe1, pe2, l1_.value(), l2_.value());

    std::array<double, NP> z{}, lumi{};
    for (int i = 0; i < NP; ++i) {
        z[i] = static_cast<double>(i + 1) / NP;
        lumi[i] =
            gagatt::photonLuminosity(z[i], X, pe1, pc1, pe2, pc2, {l1}, {l2});

        fout << z[i] << "  " << lumi[i] << '\n';
    }

    std::cout << "luminosity: the output has been stored in " << argv[1]
              << '\n';
    fout.close();
}
