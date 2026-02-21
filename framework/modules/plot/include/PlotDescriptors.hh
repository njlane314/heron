/* -- C++ -- */
/**
 *  @file  framework/plot/include/PlotDescriptors.hh
 *
 *  @brief Plot configuration descriptors that define plot metadata, labels,
 *         and binning settings for output visualisations.
 */

#ifndef HERON_PLOT_DESCRIPTORS_H
#define HERON_PLOT_DESCRIPTORS_H

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
    // STV/CCNp-style presentation: top legend band, thicker axes/ticks,
    // dotted uncertainty outline, and Events/bin defaults.
    bool stv_style = true;
    // MicroBooNE-like defaults: ratio panel + uncertainty band in ratio.
    bool show_ratio = true;
    bool show_ratio_band = true;
    // Draw counts as a density (events / bin width), matching common publication style.
    bool normalise_by_bin_width = true;
    // Draw a small χ² box on the main pad (χ² / Nbins in visible range).
    bool show_chi2 = true;
    bool use_log_x = false;
    bool use_log_y = false;
    bool annotate_numbers = true;
    bool overlay_signal = false;
    bool legend_on_top = true;
    bool show_legend = true;
    bool show_watermark = true;
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


#endif // HERON_PLOT_DESCRIPTORS_H
