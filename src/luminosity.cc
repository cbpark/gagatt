#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include "helicity.h"
#include "photon.h"

constexpr double X = 4.8;
constexpr int NP = 250;

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
    if (argc != 6) {
        std::cerr
            << "usage: ./bin/luminosity <output.dat> <pe1> <pe2> <l1> <l2>\n"
            << "  <output.dat>: output file name.\n"
            << "  <pe1> <pe2>: pols of initial electrons.\n"
            << "  <l1> <l2>: pols of colliding photons.\n";
        return EXIT_FAILURE;
    }

    auto pe1_ = parse_double(argv[2], "pe1");
    auto pe2_ = parse_double(argv[3], "pe2");
    auto l1_ = parse_double(argv[4], "l1");
    auto l2_ = parse_double(argv[5], "l2");
    if (!pe1_ || !pe2_ || !l1_ || !l2_) { return EXIT_FAILURE; }

    const double pe1 = *pe1_;
    const double pe2 = *pe2_;
    const double pc1 = -pe1;
    const double pc2 = -pe2;

    const auto l1 = gagatt::toHelicity(*l1_);
    const auto l2 = gagatt::toHelicity(*l2_);

    std::cout << "luminosity: (pe1, pe2) = (" << pe1 << ", " << pe2 << "), "
              << "(l1, l2) = (" << *l1_ << ", " << *l2_ << ")\n";

    std::ofstream fout(argv[1]);
    if (!fout) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    for (int i = 0; i < NP; ++i) {
        const double z = static_cast<double>(i + 1) / NP;
        const double lumi =
            gagatt::photonLuminosity(z, X, pe1, pc1, pe2, pc2, l1, l2);

        fout << z << "  " << lumi << '\n';
    }

    std::cout << "luminosity: the output has been stored in " << argv[1]
              << '\n';
}
