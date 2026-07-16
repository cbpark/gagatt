#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "mc.h"

using namespace gagatt;

void printUsage(const char *prog) {
    std::cerr
        << "usage: " << prog
        << " <sqrt_s> <pe1> <pe2>"
           " [n_events] [L_scan_min_ab] [L_scan_max_ab] [n_L_points]"
           " [--coscut <val>] [--sqrts_max <val>]\n\n"
        << "  sqrt_s        : e+e- CM energy [GeV] (e.g. 500 or 1000)\n"
        << "  pe1           : electron polarization (e.g. +1.0 or -1.0)\n"
        << "  pe2           : positron polarization (e.g. +1.0 or -1.0)\n"
        << "  n_events      : number of MC events [optional, default "
           "1000000]\n"
        << "  L_scan_min_ab : scan start [ab^-1] [optional, default "
           "0.001]\n"
        << "  L_scan_max_ab : scan end   [ab^-1] [optional, default "
           "1.0]\n"
        << "  n_L_points    : number of luminosity pts [optional, default "
           "200]\n"
        << "  --coscut <val>    : restrict to |cos(Theta)| < val,"
           " range (0, 1)\n"
        << "  --sqrts_max <val> : upper edge of sqrts phase-space scan"
           " [GeV]\n"
        << "\n"
        << "Benchmark scenarios:\n"
        << "  Section 5.4: sqrt_s=500  pe1=+1 pe2=+1\n"
        << "  Section 5.5: sqrt_s=1000 pe1=+1 pe2=-1\n"
        << "  Unpolarized: sqrt_s=500  pe1=0  pe2=0\n";
}

int main(int argc, char *argv[]) {
    // Strip --coscut and --sqrts_max from argv first, before positional
    // parsing.  This avoids any getopt confusion with negative numbers
    // such as pe2 = -1.
    MCConfig cfg;

    std::vector<const char *> pos;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);

        auto nextVal = [&](const char *flag) -> const char * {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: " << flag << " requires a value\n";
                std::exit(EXIT_FAILURE);
            }
            return argv[++i];
        };

        if (arg == "--coscut") {
            const double v = std::stod(nextVal("--coscut"));
            if (v <= 0.0 || v >= 1.0) {
                std::cerr << std::format(
                    "ERROR: --coscut must be in (0, 1); got {:.4f}\n", v);
                return EXIT_FAILURE;
            }
            cfg.cos_th_min = -v;
            cfg.cos_th_max = v;
        } else if (arg == "--sqrts_max") {
            cfg.sqrts_max = std::stod(nextVal("--sqrts_max"));
        } else {
            pos.push_back(argv[i]);
        }
    }

    if (pos.size() < 3 || pos.size() > 7) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    cfg.sqrt_s = std::stod(pos[0]);
    cfg.pe1 = std::stod(pos[1]);
    cfg.pe2 = std::stod(pos[2]);
    if (pos.size() >= 4) { cfg.n_events = std::stoll(pos[3]); }
    if (pos.size() >= 5) { cfg.L_scan_min_ab = std::stod(pos[4]); }
    if (pos.size() >= 6) { cfg.L_scan_max_ab = std::stod(pos[5]); }
    if (pos.size() >= 7) { cfg.n_L_points = std::stoi(pos[6]); }

    // sqrts_max defaults to sqrt_s unless --sqrts_max was given.
    if (cfg.sqrts_max < 0.0) { cfg.sqrts_max = cfg.sqrt_s; }
    if (cfg.sqrts_max > cfg.sqrt_s) {
        std::cerr << std::format(
            "ERROR: --sqrts_max ({:.1f}) cannot exceed sqrt_s ({:.1f})\n",
            cfg.sqrts_max, cfg.sqrt_s);
        return EXIT_FAILURE;
    }

    const MCResult res = runMC(cfg);

    // Write luminosity scan to a .dat file
    if (!res.lumi_scan.empty()) {
        std::string fname =
            std::format("significance_sqrts{:.0f}_pe1{:+.0f}_pe2{:+.0f}",
                        cfg.sqrt_s, cfg.pe1, cfg.pe2);
        if (cfg.cos_th_min != -1.0 || cfg.cos_th_max != 1.0) {
            fname += std::format("_coscut{:.2f}", cfg.cos_th_max);
        }
        if (cfg.sqrts_max != cfg.sqrt_s) {
            fname += std::format("_sqrts_max{:.0f}", cfg.sqrts_max);
        }
        fname += ".dat";
        std::ofstream ofs(fname);
        if (!ofs) {
            std::cerr << "ERROR: cannot open " << fname << " for writing\n";
            return EXIT_FAILURE;
        }

        ofs << std::format("# significance vs luminosity\n");
        ofs << std::format("# sqrt_s = {:.1f} GeV\n", cfg.sqrt_s);
        ofs << std::format("# pe1    = {:+.2f}\n", cfg.pe1);
        ofs << std::format("# pe2    = {:+.2f}\n", cfg.pe2);
        // Optional cuts
        if (cfg.cos_th_min != -1.0 || cfg.cos_th_max != 1.0) {
            ofs << std::format("# cos_th_min = {:.4f}\n", cfg.cos_th_min);
            ofs << std::format("# cos_th_max = {:.4f}\n", cfg.cos_th_max);
        }
        if (cfg.sqrts_max != cfg.sqrt_s) {
            ofs << std::format("# sqrts_max  = {:.1f} GeV\n", cfg.sqrts_max);
        }
        ofs << std::format("# N_MC   = {}\n", res.n_events_generated);
        ofs << std::format("# sigma_tot      = {:.6f} fb\n", res.total_xsec_fb);
        ofs << std::format("# theory C  = {:.6f}\n", res.theory_concurrence);
        ofs << std::format("# MC C      = {:.6f}\n", res.mc_concurrence);
        ofs << std::format("# theory D  = {:.6f}\n", res.theory_D);
        ofs << std::format("# MC D      = {:.6f}\n", res.mc_D);
        ofs << std::format("# theory m12= {:.6f}\n", res.theory_m12);
        ofs << std::format("# MC m12    = {:.6f}\n", res.mc_m12);
        ofs << std::format(
            "# note: significance(L) = significance(N_MC)"
            " * sqrt(sigma_tot[fb] * L[ab^-1] * 1e3 / N_MC)\n");
        ofs << std::format("#\n");
        ofs << std::format("# {:>12s} {:>16s} {:>16s} {:>16s}\n", "L[ab^-1]",
                           "sig_C[sigma]", "sig_D[sigma]", "sig_m12[sigma]");

        for (const auto &pt : res.lumi_scan) {
            ofs << std::format(" {:12.6f} {:16.6f} {:16.6f} {:16.6f}\n",
                               pt.L_ab, pt.significance_concurrence,
                               pt.significance_D, pt.significance_m12);
        }

        std::cout << std::format("\n-- scan written to: {}\n", fname);
    }

    return (res.n_events_generated > 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
