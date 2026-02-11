/* -- C++ -- */
/**
 *  @file  plot/include/PlotDescriptors.hh
 *
 *  @brief Plot configuration descriptors that define plot metadata, labels,
 *         and binning settings for output visualisations.
 */

#ifndef NUXSEC_PLOT_DESCRIPTORS_H
#define NUXSEC_PLOT_DESCRIPTORS_H

#include <cctype>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TMatrixDSym.h>

#include "../../ana/include/SelectionService.hh"


namespace nu
{

enum class CutDir
{
    LessThan,
    GreaterThan
};

struct CutSpec
{
    double x = 0.0;
    CutDir dir = CutDir::GreaterThan;
};

struct Entry
{
    SelectionEntry selection;
    double pot_nom = 0.0;
    double pot_eqv = 0.0;
    std::string beamline;
    std::string period;

    ROOT::RDF::RNode rnode() const { return selection.nominal.rnode(); }
};

struct Options
{
    std::string out_dir = ".";
    std::string image_format = "png";
    bool show_ratio = false;
    bool show_ratio_band = false;
    bool use_log_y = false;
    bool annotate_numbers = false;
    bool overlay_signal = false;
    bool legend_on_top = false;
    bool show_legend = true;
    bool show_cuts = false;
    double legend_split = 0.75;
    double y_min = 0.0;
    double y_max = -1.0;
    std::string x_title;
    std::string y_title;
    std::vector<int> signal_channels;
    std::shared_ptr<TMatrixDSym> total_cov;
    std::vector<double> syst_bin;
    std::vector<CutSpec> cuts;
    double total_protons_on_target = 0.0;
    std::string beamline;
    std::vector<std::string> run_numbers;
    std::vector<std::string> periods;
    std::string analysis_region_label;
    // --- Adaptive ("minimum-stat-per-bin") binning ---
    bool adaptive_binning = false;
    double adaptive_min_sumw = 0.0;       // Wmin (sum of weights; <=0 disables)
    double adaptive_max_relerr = 0.0;     // relErrMax (sqrt(sumw2)/|sumw|; <=0 disables)
    bool adaptive_fold_overflow = true;   // fold under/overflow before edge-making + rebin

    // Keep this many *fine* bins fixed (unmerged) at each edge of the plot.
    // This gives you constant-width leading/trailing bins (often intentionally empty padding).
    int adaptive_edge_bins = 0;

    // IMPORTANT: the min-stat algorithm only MERGES bins. If you fill with your final coarse
    // binning, results will look uniform. Fill finer, then merge.
    // Filled nbins = spec.nbins * adaptive_fine_bin_factor (clamped).
    int adaptive_fine_bin_factor = 10;

    // --- Particle-level plotting (vector-valued branches) ---
    // If enabled, MC is stacked by truth-matched particle type using `particle_pdg_branch`.
    // This expects `spec.expr` (or `spec.id` if expr is empty) to evaluate to a vector-like
    // column (e.g. track_length, backtracked_energies, ...). Data is still drawn unstacked.
    bool particle_level = false;
    std::string particle_pdg_branch = "backtracked_pdg_codes";
    bool particle_drop_nan = true;

    // Channel column used for stack/unstack grouping (default analysis channel branch).
    std::string channel_column = "analysis_channels";

    // Optional channel overrides for specialised overlays (e.g. occupancy components).
    // If empty, built-in Channels defaults are used.
    std::vector<int> unstack_channel_keys;
    std::map<int, std::string> unstack_channel_labels;
    std::map<int, int> unstack_channel_colours;
};

struct TH1DModel
{
    std::string id;
    std::string name;
    std::string title;
    std::string expr;
    std::string weight = "w_nominal";
    int nbins = 1;
    double xmin = 0.0;
    double xmax = 1.0;
    Preset sel = Preset::Muon;

    ROOT::RDF::TH1DModel model(const std::string &suffix = "") const
    {
        const std::string base = !id.empty() ? id : name;
        const std::string hist_name = sanitise(base + suffix);
        const std::string hist_title = title.empty() ? base : title;
        return ROOT::RDF::TH1DModel(hist_name.c_str(), hist_title.c_str(), nbins, xmin, xmax);
    }

    std::string axis_title() const
    {
        if (!title.empty())
        {
            return title;
        }
        const std::string base = !name.empty() ? name : id;
        if (base.empty())
        {
            return ";x;Events";
        }
        return ";" + base + ";Events";
    }

  private:
    static std::string sanitise(const std::string &raw)
    {
        std::string out;
        out.reserve(raw.size());
        for (unsigned char c : raw)
        {
            if (std::isalnum(c) || c == '_' || c == '-')
            {
                out.push_back(static_cast<char>(c));
            }
            else if (c == '.' || c == ' ')
            {
                out.push_back('_');
            }
            else
            {
                out.push_back('_');
            }
        }
        if (out.empty())
        {
            return "plot";
        }
        return out;
    }
};

} // namespace nu


#endif // NUXSEC_PLOT_DESCRIPTORS_H
