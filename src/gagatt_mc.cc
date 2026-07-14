#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include "mc.h"

using namespace gagatt;

void printUsage(const char *prog) {
    std::cerr
        << "usage: " << prog
        << " <sqrt_s> <pe1> <pe2> <n_events>"
           " [L_scan_min_ab] [L_scan_max_ab] [n_L_points]\n"
        << "\n"
        << "  sqrt_s        : e+e- CM energy [GeV]      (e.g. 500 or 1000)\n"
        << "  pe1           : electron polarization     (e.g. +1.0 or -1.0)\n"
        << "  pe2           : positron polarization     (e.g. +1.0 or -1.0)\n"
        << "  n_events      : number of MC events       (e.g. 1000000)\n"
        << "  L_scan_min_ab : scan start [ab^-1]        [optional, default "
           "0.01]\n"
        << "  L_scan_max_ab : scan end   [ab^-1]        [optional, default "
           "2.0]\n"
        << "  n_L_points    : number of luminosity pts  [optional, default "
           "200]\n"
        << "\n"
        << "Benchmark scenarios:\n"
        << "  Section 5.4: sqrt_s=500  pe1=+1 pe2=+1\n"
        << "  Section 5.5: sqrt_s=1000 pe1=+1 pe2=-1\n"
        << "  Unpolarized: sqrt_s=500  pe1=0  pe2=0\n";
}

int main(int argc, char *argv[]) {
    if (argc < 5 || argc > 8) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    MCConfig cfg;
    cfg.sqrt_s = std::stod(argv[1]);
    cfg.pe1 = std::stod(argv[2]);
    cfg.pe2 = std::stod(argv[3]);
    cfg.n_events = std::stoll(argv[4]);

    // Luminosity scan parameters (optional)
    cfg.L_scan_min_ab = (argc >= 6) ? std::stod(argv[5]) : 0.01;
    cfg.L_scan_max_ab = (argc >= 7) ? std::stod(argv[6]) : 2.0;
    cfg.n_L_points = (argc >= 8) ? std::stoi(argv[7]) : 200;

    // L_ee_fb is only used internally by mc.cc for the weight normalization;
    // set it to n_events so weights come out in natural units.
    // The significance scan uses sigma_tot and N_MC directly:see mc.cc Phase 7.
    cfg.L_ee_fb = static_cast<double>(cfg.n_events);  // dummy; see mc.cc
    cfg.sqrts_max = cfg.sqrt_s;

    const MCResult res = runMC(cfg);

    // ------------------------------------------------------------------
    // Write luminosity scan to a .dat file
    // ------------------------------------------------------------------
    if (!res.lumi_scan.empty()) {
        const std::string fname =
            std::format("significance_sqrts{:.0f}_pe1{:+.0f}_pe2{:+.0f}.dat",
                        cfg.sqrt_s, cfg.pe1, cfg.pe2);

        std::ofstream ofs(fname);
        if (!ofs) {
            std::cerr << "ERROR: cannot open " << fname << " for writing\n";
            return EXIT_FAILURE;
        }

        // Header
        ofs << std::format("# significance vs luminosity\n");
        ofs << std::format("# sqrt_s    = {:.1f} GeV\n", cfg.sqrt_s);
        ofs << std::format("# pe1       = {:+.2f}\n", cfg.pe1);
        ofs << std::format("# pe2       = {:+.2f}\n", cfg.pe2);
        ofs << std::format("# N_MC      = {}\n", res.n_events_generated);
        ofs << std::format("# sigma_tot = {:.6f} fb\n", res.total_xsec_fb);
        ofs << std::format("# theory m12= {:.6f}\n", res.theory_m12);
        ofs << std::format("# MC m12    = {:.6f}\n", res.mc_m12);
        ofs << std::format(
            "# note: significance(L) = significance(N_MC)"
            " * sqrt(sigma_tot[fb] * L[fb^-1] / N_MC)\n");
        ofs << std::format("#\n");
        ofs << std::format("# {:>12s}  {:>16s}  {:>16s}\n", "L[ab^-1]",
                           "sig_D[sigma]", "sig_Bell[sigma]");

        for (const auto &pt : res.lumi_scan)
            ofs << std::format("  {:12.6f}  {:16.6f}  {:16.6f}\n", pt.L_ab,
                               pt.significance_D, pt.significance_bell);

        std::cout << std::format("\n-- scan written to: {}\n", fname);
    }

    return (res.n_events_generated > 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
