/* -- C++ -- */
/**
 *  @file  plot/include/AdaptiveBinningService.hh
 *
 *  @brief Service for deriving adaptive min-stat bin edges and rebinning histograms.
 */

#ifndef HERON_PLOT_ADAPTIVE_BINNING_SERVICE_H
#define HERON_PLOT_ADAPTIVE_BINNING_SERVICE_H

#include <cstddef>
#include <memory>
#include <string_view>
#include <vector>

class TH1D;

namespace nu
{

struct Options;

class AdaptiveBinningService
{
  public:
    struct MinStatConfig
    {
        bool enabled = false;

        /// Bin validity: require |sumw| >= min_sumw (<=0 disables).
        double min_sumw = 0.0;

        /// Bin validity: require sqrt(sumw2)/|sumw| <= max_rel_err (<=0 disables).
        double max_rel_err = 0.0;

        /// Fold under/overflow into first/last in-range bins before edge-making/rebinning.
        bool fold_overflow = true;

        /// Use |sumw| in checks (recommended when negative weights are possible).
        bool use_abs_sumw = true;

        // Keep N *fine* bins fixed on each edge (unmerged).
        int edge_bins = 0;
    };

    static AdaptiveBinningService &instance();

    static MinStatConfig config_from(const Options &opt);

    std::vector<double> edges_min_stat(const TH1D &fine, const MinStatConfig &cfg) const;

    // Preferred API: use cfg (so edge_bins + overflow folding behave consistently).
    std::unique_ptr<TH1D> rebin_to_edges(const TH1D &h,
                                         const std::vector<double> &edges,
                                         std::string_view new_name,
                                         const MinStatConfig &cfg) const;

    // Backwards-compatible wrapper.
    std::unique_ptr<TH1D> rebin_to_edges(const TH1D &h,
                                         const std::vector<double> &edges,
                                         std::string_view new_name,
                                         bool fold_overflow) const;

    // Preferred API: use cfg (so edge_bins + overflow folding behave consistently).
    std::unique_ptr<TH1D> sum_hists(const std::vector<const TH1D *> &parts,
                                    std::string_view new_name,
                                    const MinStatConfig &cfg) const;

    // Backwards-compatible wrapper.
    std::unique_ptr<TH1D> sum_hists(const std::vector<const TH1D *> &parts,
                                    std::string_view new_name,
                                    bool fold_overflow) const;

    void fold_overflow(TH1D &h) const;

  private:
    AdaptiveBinningService() = default;
};

} // namespace nu

#endif // HERON_PLOT_ADAPTIVE_BINNING_SERVICE_H
