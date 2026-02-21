/* -- C++ -- */
/**
 *  @file  framework/ana/include/TemplateBinningOptimizer1D.hh
 *
 *  @brief Greedy adjacent-bin merger for 1D reco-space binning, optimised for
 *         forward-folded binned likelihood/template fits.
 */

#ifndef HERON_ANA_TEMPLATE_BINNING_OPTIMIZER_1D_H
#define HERON_ANA_TEMPLATE_BINNING_OPTIMIZER_1D_H

#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

class TH1;

/**
 *  @brief Template-based 1D binning optimiser.
 *
 *  Core ideas:
 *   - Objective: expected POI precision (Fisher / profiled Fisher)
 *   - Constraints: minimum expected yield per bin and/or maximum relative MC
 *                  statistical uncertainty per bin
 *
 *  Inputs are ROOT histograms (nominal + parameter variations or derivatives).
 *  Output is a vector of bin edges suitable for TH1::Rebin(..., edges.data()).
 */
class TemplateBinningOptimizer1D final
{
  public:
    struct Parameter
    {
        std::string name;

        /// Provide either a derivative histogram, or up/down templates + step.
        const TH1 *p_derivative = nullptr;
        const TH1 *p_up = nullptr;
        const TH1 *p_down = nullptr;
        double step = 1.0;

        /// Gaussian prior width for nuisance parameters (0 => unconstrained/free).
        double prior_sigma = 0.0;

        /// Exactly one parameter should be the POI for the binning objective.
        bool is_poi = false;
    };

    struct Channel
    {
        std::string name;
        /// Expected total (signal+background) reco template, fine-binned.
        const TH1 *p_nominal = nullptr;
        /// POI + selected nuisances (optional).
        std::vector<Parameter> parameters;
    };

    struct Config
    {
        /// Constraints applied during merging.
        double mu_min = 1.0;
        double rel_mc_max = 0.15;

        /// If true, every channel must satisfy constraints in each bin.
        bool require_per_channel_constraints = true;

        /// Objective behaviour.
        bool profile_nuisances = true;
        double mu_floor_for_objective = 1e-12;

        /// Keep merging until number of bins <= max_bins if max_bins > 0.
        int max_bins = -1;

        /// Optional small preference against wide bins (tie-breaker).
        double width_penalty = 0.0;

        /// Diagnostics.
        bool verbose = false;
        std::ostream *p_log = nullptr;
    };

    struct BinReport
    {
        double low = 0.0;
        double high = 0.0;

        double mu_sum = 0.0;
        double rel_mc_worst = 0.0;

        bool passes_constraints = true;
    };

    struct Result
    {
        std::vector<double> edges;
        double expected_sigma_poi = std::numeric_limits<double>::infinity();
        std::vector<BinReport> bins;
    };

    explicit TemplateBinningOptimizer1D(Config cfg);

    /// Optimise using a single channel.
    Result optimise(const Channel &channel) const;

    /// Optimise using multiple channels and return shared edges.
    Result optimise(const std::vector<Channel> &channels) const;

  private:
    Config m_cfg;
};

#endif // HERON_ANA_TEMPLATE_BINNING_OPTIMIZER_1D_H
