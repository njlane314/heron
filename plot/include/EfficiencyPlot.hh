/* -- C++ -- */
/**
 *  @file  plot/include/EfficiencyPlot.hh
 *
 *  @brief Efficiency plotting helper that builds denominator/numerator
 *         distributions and overlays right-axis efficiency.
 */

#ifndef HERON_PLOT_EFFICIENCY_PLOT_H
#define HERON_PLOT_EFFICIENCY_PLOT_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>

#include "Plotter.hh"
#include "PlottingHelper.hh"


namespace nu
{

class EfficiencyPlot
{
  public:
    struct Config
    {
        std::string x_title;
        std::string y_counts_title = "Events";
        std::string y_eff_title = "Efficiency";

        std::string legend_total = "All events";
        std::string legend_passed = "Selected";
        std::string legend_eff = "Efficiency";

        bool draw_distributions = true;
        bool draw_total_hist = true;
        bool draw_passed_hist = true;
        bool logy = false;

        double eff_ymin = 0.0;
        double eff_ymax = 1.05;

        bool auto_x_range = false;
        double x_pad_fraction = 0.05;

        double conf_level = 0.68;
        TEfficiency::EStatOption stat = TEfficiency::kFCP;
        bool use_weighted_events = false;

        // --- Presentation / logging ---
        // Do NOT draw denom/numer/overall eff text on the plot (cleaner, publication-like).
        // Instead print to stdout (see print_stats).
        bool draw_stats_text = false;

        // Print denom/numer/overall eff for each saved plot.
        bool print_stats = true;

        // Allow scientific-notation "×10^{n}" on the left (counts) axis by default.
        bool no_exponent_y = false;

        std::vector<std::string> extra_text_lines;
    };

    EfficiencyPlot(TH1DModel spec, Options opt = Options{});
    EfficiencyPlot(TH1DModel spec, Options opt, Config cfg);

    const TH1DModel &spec() const noexcept { return spec_; }
    const Options &options() const noexcept { return opt_; }

    Config &config() noexcept { return cfg_; }
    const Config &config() const noexcept { return cfg_; }

    int compute(ROOT::RDF::RNode denom_node, ROOT::RDF::RNode pass_node);

    int compute(ROOT::RDF::RNode base,
                const std::string &denom_sel,
                const std::string &pass_sel,
                const std::string &extra_sel = "true");

    int draw_and_save(const std::string &file_stem,
                      const std::string &format = "") const;

    const TH1D *total_hist() const noexcept { return h_total_.get(); }
    const TH1D *passed_hist() const noexcept { return h_passed_.get(); }
    const TGraphAsymmErrors *eff_graph() const noexcept { return g_eff_.get(); }

    std::uint64_t denom_entries() const noexcept { return n_denom_; }
    std::uint64_t pass_entries() const noexcept { return n_pass_; }

    bool ready() const noexcept { return ready_; }

  private:
    TH1DModel spec_;
    Options opt_;
    Config cfg_;

    std::unique_ptr<TH1D> h_total_;
    std::unique_ptr<TH1D> h_passed_;
    std::unique_ptr<TGraphAsymmErrors> g_eff_;

    std::uint64_t n_denom_ = 0;
    std::uint64_t n_pass_ = 0;

    bool ready_ = false;

    static std::string sanitise_(const std::string &s);
};

} // namespace nu


#endif // HERON_PLOT_EFFICIENCY_PLOT_H
