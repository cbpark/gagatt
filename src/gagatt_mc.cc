#include <cstdlib>
#include <iostream>
#include "mc.h"

using namespace gagatt;

static void printUsage(const char *prog) {
    std::cerr << "usage: " << prog
              << " <sqrt_s> <pe1> <pe2> <L_ee_fb> [n_events]\n"
              << "  sqrt_s    : e+e- CM energy [GeV]  (e.g. 500 or 1000)\n"
              << "  pe1       : electron polarization  (e.g. +1.0 or -1.0)\n"
              << "  pe2       : positron polarization  (e.g. +1.0 or -1.0)\n"
              << "  L_ee_fb   : integrated luminosity [fb^-1]  (e.g. 1000)\n"
              << "  n_events  : MC events [optional, default 1000000]\n\n"
              << "Benchmark scenarios from the paper:\n"
              << "  Section 5.4: sqrt_s=500  pe1=+1 pe2=+1\n"
              << "  Section 5.5: sqrt_s=1000 pe1=+1 pe2=-1\n"
              << "  Unpolarized: sqrt_s=500  pe1=0  pe2=0\n";
}

int main(int argc, char *argv[]) {
    if (argc < 5 || argc > 6) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    MCConfig cfg;
    cfg.sqrt_s = std::stod(argv[1]);
    cfg.pe1 = std::stod(argv[2]);
    cfg.pe2 = std::stod(argv[3]);
    cfg.L_ee_fb = std::stod(argv[4]);
    if (argc == 6) { cfg.n_events = std::stoll(argv[5]); }

    // Derived: sqrts_max = sqrt_s (kinematic ceiling)
    cfg.sqrts_max = cfg.sqrt_s;

    // Print run summary
    const MCResult res = runMC(cfg);
}
