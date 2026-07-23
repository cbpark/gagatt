#ifndef SRC_MC_HELPER_H
#define SRC_MC_HELPER_H

#include <random>
#include <utility>
#include <vector>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <Eigen/Dense>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include "mc.h"
#include "photon.h"
#include "spin_density.h"

namespace gagatt {
// -----------------------------------------------------------------------
// ZCacheEntry / LumiCache
//
// Precomputed (LumiWeights, L_tot) per sqrt_s_hat bin, indexed by j in
// [0, n_sqrts). Reused by buildWeightTable() so lumiWeightsAndTotal() is
// evaluated exactly once per sqrt_s_hat bin instead of once per (i, j) cell.
// -----------------------------------------------------------------------
struct ZCacheEntry {
    LumiWeights lw;
    double L_tot;
};

// Builds the z-cache over n_sqrts bins spanning [sqrts_min, sqrts_max].
// pc1 = -pe1, pc2 = -pe2.
std::vector<ZCacheEntry> buildLumiCache(const MCConfig &cfg, double pc1,
                                        double pc2, double sqrts_min,
                                        double d_sqrts);

// -----------------------------------------------------------------------
// partialXsec:
//
// differential partonic cross section d sigma_hat / d cos_th [GeV^-2]
// = (beta Nc / 32 pi s_hat) * |A_C|^2 * (C_1^w + C_3^w)
//
// sdc.norm_factor = C_1^w + C_3^w, with |A_C|^2 already absorbed via
// overall_fac^2 = (COUPLING_FACTOR / denom)^2 inside polCoeffsForHelicity.
// For unpolarized photon average, divide by 4 afterwards.
// This was previously static in mc_helper.cc – now public so that
// sqrt_s_hat_xsec.cc and other analyses can reuse it.
// -----------------------------------------------------------------------
double partialXsec(double sqrt_s_hat, const SDMatrixCoefficients &sdc);

// -----------------------------------------------------------------------
// eventRate
//
// differential event rate d^2 sigma / (d sqrt_s_hat d cos_th) [fb/GeV].
// See mc_helper.cc for the full derivation/unit bookkeeping.
// -----------------------------------------------------------------------
double eventRate(double sqrt_s_hat, const SDMatrixCoefficients &sdc,
                 double L_tot, double sqrt_s);

// -----------------------------------------------------------------------
// WeightTable
//
// Flattened 2-D (cos_th, sqrt_s_hat) weight table plus the cached
// SDMatrixCoefficients for each bin (so the event loop doesn't need to
// recompute sdc for the sampled bin), and the luminosity+phase-space
// weighted theory averages of concurrence, D,  m12.
// -----------------------------------------------------------------------
struct WeightTable {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    std::vector<double> bin_weights;
    std::vector<SDMatrixCoefficients> sdc_cache;
    double total_weight = 0.0;

    double theory_concurrence = 0.0;
    double theory_D = 0.0;
    double theory_m12 = 0.0;
};

// Builds the flattened weight table over the (n_cos x n_sqrts) grid and
// accumulates the theory averages in a single pass. Returns an empty
// bin_weights vector (size 0) if total_weight <= 0 (caller must check).
WeightTable buildWeightTable(const MCConfig &cfg,
                             const std::vector<ZCacheEntry> &zcache,
                             double sqrts_min, double d_sqrts, double d_cos);

// -----------------------------------------------------------------------
// sampleDecayAngles
//
// Draw (q+, q-) unit vectors for l+/l- in the top/anti-top rest frames
// from the joint angular distribution via accept/reject sampling.
// See mc_helper.cc for the distribution and envelope details.
// -----------------------------------------------------------------------
std::pair<Eigen::Vector3d, Eigen::Vector3d> sampleDecayAngles(
    const SDMatrixCoefficients &sdc, std::mt19937_64 &rng);

// Per-bin moment accumulator (one entry per (cos_th, sqrt_s_hat) bin).
struct BinMoments {
    Eigen::Matrix3d S1_qpqm = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d S2_qpqm = Eigen::Matrix3d::Zero();
    Eigen::Vector3d S1_qp = Eigen::Vector3d::Zero();
    Eigen::Vector3d S2_qp = Eigen::Vector3d::Zero();
    Eigen::Vector3d S1_qm = Eigen::Vector3d::Zero();
    Eigen::Vector3d S2_qm = Eigen::Vector3d::Zero();
    long long n = 0;
};

// -----------------------------------------------------------------------
// EventLoopResult
//
// Accumulated event count and per-bin moments from the MC event loop.
// All sigma estimation is done from per-bin S2_qpqm in BinMoments.
// -----------------------------------------------------------------------
struct EventLoopResult {
    long long n_accepted = 0;

    std::vector<BinMoments> per_bin;
};

// Runs the MC event loop: draws a (cos_th, sqrt_s_hat) bin proportional to
// wt.bin_weights, samples decay angles, and accumulates moments.
EventLoopResult runEventLoop(long long n_events, const WeightTable &wt,
                             std::mt19937_64 &rng, bool verbose = false);

// -----------------------------------------------------------------------
// ReconstructedMC
//
// All MC-side quantities derived from the raw event-loop moments:
// B_+, B_-, C_ij, their uncertainties, D, concurrence, and the
// Horodecki m12 parameter.
// -----------------------------------------------------------------------
struct ReconstructedMC {
    double mc_concurrence = 0.0;
    double sigma_concurrence = 0.0;
    double significance_concurrence = 0.0;

    double mc_D = 0.0;
    double sigma_D = 0.0;
    double significance_D = 0.0;  // (-1/3 - D) / sigma_D when D < -1/3

    double mc_m12 = 0.0;
    double sigma_m12 = 0.0;
    double significance_m12 = 0.0;  // (m12 - 1) / sigma_m12 when m12 > 1
};

// Turns raw event-loop moments into C_ij and all derived quantities.
ReconstructedMC reconstructFromMoments(const EventLoopResult &ev);

// -----------------------------------------------------------------------
// Luminosity scan
//
// At each luminosity point L, compute N(L) = sigma_eff_fb * L_fb, then
// run n_lumi_seeds independent MC event loops each with N(L) events,
// reconstruct the observables from each seed, and report the mean and
// std-dev of the significance across seeds.
//
// This correctly captures the degradation of mc_concurrence itself at
// low luminosity (bin-emptying effect), which the old sqrt(N) rescaling
// shortcut could not account for.
// -----------------------------------------------------------------------
std::vector<LumiScanPoint> computeLumiScan(const MCConfig &cfg,
                                           double sigma_eff_fb,
                                           const WeightTable &wt,
                                           std::mt19937_64 &rng);

// Integrate over cos_th in [cos_min, cos_max] using Simpson's rule.
template <typename Func>
double integrate_cos(Func f, double cos_min, double cos_max, int n) {
    if (n % 2 == 1) ++n;
    const double d = (cos_max - cos_min) / n;
    double sum = f(cos_min) + f(cos_max);
    for (int i = 1; i < n; ++i) {
        const double c = cos_min + i * d;
        const double fi = f(c);
        sum += (i % 2 == 0) ? 2.0 * fi : 4.0 * fi;
    }
    return sum * d / 3.0;
}
}  // namespace gagatt

#endif  // SRC_MC_HELPER_H
