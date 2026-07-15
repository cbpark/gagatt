#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include "mc.h"
#include "mc_helper.h"

int main(int argc, char *argv[]) {
    // ------------------------------------------------------------------
    // Usage: gagatt_dscan <sqrt_s> <pe1> <pe2> [options]
    //
    //   <sqrt_s>          CM energy [GeV], e.g. 500
    //   <pe1>             electron beam polarization
    //   <pe2>             positron beam polarization
    //
    // Optional (with defaults):
    //   --bins <n>        number of sqrt(s_hat) scan bins (default: 50)
    //   --ncos <n>        cos(Theta) sub-bins per scan bin (default: 50)
    //   --nevents <n>     MC events per bin (default: 100000)
    //   --cosmin <v>      cos(Theta) cut lower bound (default: -1.0)
    //   --cosmax <v>      cos(Theta) cut upper bound (default:  1.0)
    //   --seed <n>        RNG seed (default: 0 = random_device)
    //   --output <file>   output .dat file (default: dscan.dat)
    // ------------------------------------------------------------------
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <sqrt_s> <pe1> <pe2> [options]\n";
        return EXIT_FAILURE;
    }

    // Positional arguments
    const double sqrt_s = std::stod(argv[1]);
    const double pe1 = std::stod(argv[2]);
    const double pe2 = std::stod(argv[3]);

    // Defaults
    int n_scan_bins = 50;
    int n_cos_sub = 50;
    long long n_events_per_bin = 100000;
    double cos_th_min = -1.0;
    double cos_th_max = 1.0;
    uint64_t seed = 0;
    std::string outfile = "dscan.dat";

    // Parse optional keyword arguments
    for (int i = 4; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--bins" && i + 1 < argc) {
            n_scan_bins = std::stoi(argv[++i]);
        } else if (arg == "--ncos" && i + 1 < argc) {
            n_cos_sub = std::stoi(argv[++i]);
        } else if (arg == "--nevents" && i + 1 < argc) {
            n_events_per_bin = std::stoll(argv[++i]);
        } else if (arg == "--cosmin" && i + 1 < argc) {
            cos_th_min = std::stod(argv[++i]);
        } else if (arg == "--cosmax" && i + 1 < argc) {
            cos_th_max = std::stod(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            seed = std::stoull(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            outfile = argv[++i];
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return EXIT_FAILURE;
        }
    }

    // Validate
    if (cos_th_min >= cos_th_max) {
        std::cerr << "Error: cosmin must be < cosmax.\n";
        return EXIT_FAILURE;
    }
    if (n_scan_bins <= 0 || n_cos_sub <= 0 || n_events_per_bin <= 0) {
        std::cerr << "Error: bins, ncos, and nevents must be positive.\n";
        return EXIT_FAILURE;
    }

    // Build MCConfig — sqrts_min/max are the physical threshold / sqrt_s
    gagatt::MCConfig cfg;
    cfg.sqrt_s = sqrt_s;
    cfg.pe1 = pe1;
    cfg.pe2 = pe2;
    cfg.cos_th_min = cos_th_min;  // used inside runDScanVsSqrtS
    cfg.cos_th_max = cos_th_max;
    cfg.seed = static_cast<int>(seed);
    // sqrts_min / sqrts_max are read from MCConfig defaults (2*M_top / sqrt_s)
    // — confirm these are set correctly in MCConfig's constructor/defaults.

    // Build DScanConfig
    gagatt::DScanConfig scan_cfg;
    scan_cfg.n_scan_bins = n_scan_bins;
    scan_cfg.n_cos_sub = n_cos_sub;
    scan_cfg.n_events_per_bin = n_events_per_bin;
    scan_cfg.cos_th_min = cos_th_min;
    scan_cfg.cos_th_max = cos_th_max;

    const uint64_t actual_seed =
        (seed == 0) ? static_cast<uint64_t>(std::random_device{}()) : seed;

    std::cout << std::format(
        "gagatt_dscan: sqrt_s={:.0f} GeV, pe1={:+.2f}, pe2={:+.2f}\n", sqrt_s,
        pe1, pe2);
    std::cout << std::format("  {} scan bins, {} cos sub-bins, {} events/bin\n",
                             n_scan_bins, n_cos_sub, n_events_per_bin);
    std::cout << std::format("  cos(Theta) cut: [{:.2f}, {:.2f}]\n", cos_th_min,
                             cos_th_max);

    // Run the scan
    const std::vector<gagatt::DScanBinResult> results =
        gagatt::runDScanVsSqrtS(cfg, scan_cfg, actual_seed);

    // Write output .dat file
    std::ofstream ofs(outfile);
    if (!ofs) {
        std::cerr << "Error: cannot open output file: " << outfile << "\n";
        return EXIT_FAILURE;
    }

    // Header
    ofs << std::format(
        "# gagatt_dscan: sqrt_s={:.0f} GeV, "
        "pe1={:+.2f}, pe2={:+.2f}\n",
        sqrt_s, pe1, pe2);
    ofs << std::format("# cos(Theta) cut: [{:.2f}, {:.2f}]\n", cos_th_min,
                       cos_th_max);
    ofs << std::format(
        "# {} scan bins, {} cos sub-bins, "
        "{} events/bin\n",
        n_scan_bins, n_cos_sub, n_events_per_bin);
    ofs << std::format("# seed: {}\n", actual_seed);
    ofs << std::format("#\n");
    ofs << std::format(
        "# {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} "
        "{:>12s} {:>12s} {:>14s}\n",
        "sqrts_lo", "sqrts_hi", "sqrts_mid", "theory_D", "mc_D", "sigma_D",
        "sig_D", "bin_xsec_fb");

    // Data rows
    for (const auto &res : results) {
        ofs << std::format(
            "  {:12.4f} {:12.4f} {:12.4f} {:12.6f} {:12.6f} {:12.6f} "
            "{:12.4f} {:14.6e}\n",
            res.sqrts_lo, res.sqrts_hi, res.sqrts_mid, res.theory_D, res.mc_D,
            res.sigma_D, res.significance_D, res.bin_xsec_fb);
    }

    std::cout << std::format("-- wrote {} rows to {}\n", results.size(),
                             outfile);
    return EXIT_SUCCESS;
}
