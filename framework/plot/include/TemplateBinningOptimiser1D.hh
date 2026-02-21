/* -- C++ -- */
/**
 *  @file  plot/include/TemplateBinningOptimiser1D.hh
 *
 *  @brief Greedy adjacent-bin merger for 1D reco-space binning, optimised for
 *         forward-folded binned likelihood/template fits.
 */

#ifndef HERON_PLOT_TEMPLATE_BINNING_OPTIMISER1D_H
#define HERON_PLOT_TEMPLATE_BINNING_OPTIMISER1D_H

#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

#include <TH1.h>
#include <TMatrixDSym.h>


namespace nu
{
namespace binning
{

class TemplateBinningOptimiser1D
{
  public:
    struct Parameter
    {
        std::string name;

        /// Provide either a derivative histogram or up/down templates.
        const TH1 *derivative = nullptr;
        const TH1 *up = nullptr;
        const TH1 *down = nullptr;
        double step = 1.0;

        /// Gaussian prior width for nuisance parameters (0 => unconstrained/free).
        double prior_sigma = 0.0;

        /// Exactly one parameter should be the POI for the binning objective.
        bool is_poi = false;
    };

    struct Channel
    {
        std::string name;
        const TH1 *nominal = nullptr;
        std::vector<Parameter> parameters;
    };

    struct Config
    {
        double mu_min = 1.0;
        double rel_mc_max = 0.15;

        bool require_per_channel_constraints = true;

        bool profile_nuisances = true;
        double mu_floor_for_objective = 1e-12;

        int max_bins = -1;

        double width_penalty = 0.0;

        bool verbose = false;
        std::ostream *log = nullptr;
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

    explicit TemplateBinningOptimiser1D(Config cfg);

    Result optimise(const Channel &channel) const;
    Result optimise(const std::vector<Channel> &channels) const;

  private:
    Config cfg_;
};

} // namespace binning
} // namespace nu


#endif // HERON_PLOT_TEMPLATE_BINNING_OPTIMISER1D_H
