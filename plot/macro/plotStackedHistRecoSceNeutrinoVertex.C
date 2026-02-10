// plot/macro/plotStackedHistRecoSceNeutrinoVertex.C
//
// Run with:
//   ./nuxsec macro plotStackedHistRecoSceNeutrinoVertex.C
//   ./nuxsec macro plotStackedHistRecoSceNeutrinoVertex.C 'plotStackedHistRecoSceNeutrinoVertex("/path/to/samples.tsv","true",false)'
//
// Notes:
//   - This macro loads aggregated samples (samples.tsv -> SampleIO -> original analysis tree)
//   - It runs your analysis column derivations so that "analysis_channels" exists for stacking.
//   - The stack is grouped by "analysis_channels"; expr controls the x-axis variable only.
//   - MC yields are scaled by w_nominal unless an alternative weight is provided.
//   - Output dir/format follow PlotEnv defaults (NUXSEC_PLOT_DIR / NUXSEC_PLOT_FORMAT).
//   - Default input uses the generated event list (event_list_<analysis>.root).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TFile.h>

#include "AnalysisConfigService.hh"
#include "ColumnDerivationService.hh"
#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "RDataFrameService.hh"
#include "SampleCLI.hh"
#include "SampleIO.hh"


using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &p)
{
    const auto n = p.size();
    if (n < 5 || p.substr(n - 5) != ".root")
        return false;

    std::unique_ptr<TFile> f(TFile::Open(p.c_str(), "READ"));
    if (!f || f->IsZombie())
        return false;

    const bool has_refs = (f->Get("sample_refs") != nullptr);
    const bool has_events_tree = (f->Get("events") != nullptr);
    const bool has_event_tree_key = (f->Get("event_tree") != nullptr);

    return has_refs && (has_events_tree || has_event_tree_key);
}


bool debug_enabled()
{
    const char *env = std::getenv("NUXSEC_DEBUG_PLOT_STACK");
    return env != nullptr && std::string(env) != "0";
}

void debug_log(const std::string &msg)
{
    if (!debug_enabled())
    {
        return;
    }
    std::cout << "[plotStackedHistRecoSceNeutrinoVertex][debug] " << msg << "\n";
    std::cout.flush();
}
} // namespace

int plot_stacked_hist_impl(const std::string &samples_tsv,
                           const std::string &mc_weight,
                           const std::string &extra_sel,
                           bool use_logy,
                           bool include_data)
{
    ROOT::EnableImplicitMT();

    debug_log("starting plot_stacked_hist_impl");
    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotStackedHistRecoSceNeutrinoVertex] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotStackedHistRecoSceNeutrinoVertex] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    std::cout << "[plotStackedHistRecoSceNeutrinoVertex] mode=event_list\n";
    debug_log("validated input root file and entering event-list mode");

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

    debug_log("mask sizes: ext=" + std::to_string(mask_ext ? mask_ext->size() : 0) +
              ", mc=" + std::to_string(mask_mc ? mask_mc->size() : 0) +
              ", data=" + std::to_string(mask_data ? mask_data->size() : 0));

    auto filter_by_mask = [](ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask) {
        return n.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base = rdf;
    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext);
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0
                                                && sid < static_cast<int>(mask_ext->size())
                                                && (*mask_ext)[static_cast<size_t>(sid)]);
                                   },
                                   {"sample_id"});
    ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

    std::vector<Entry> entries;
    entries.reserve(include_data ? 3 : 2);

    std::vector<const Entry *> mc;
    std::vector<const Entry *> data;

    ProcessorEntry rec_mc;
    rec_mc.source = Type::kMC;

    ProcessorEntry rec_ext;
    rec_ext.source = Type::kExt;

    entries.emplace_back(make_entry(std::move(node_mc), rec_mc));
    Entry &e_mc = entries.back();
    mc.push_back(&e_mc);

    entries.emplace_back(make_entry(std::move(node_ext), rec_ext));
    Entry &e_ext = entries.back();
    mc.push_back(&e_ext);

    Entry *p_data = nullptr;
    if (include_data)
    {
        ProcessorEntry rec_data;
        rec_data.source = Type::kData;
        entries.emplace_back(make_entry(std::move(node_data), rec_data));
        p_data = &entries.back();
        data.push_back(p_data);
    }

    if (!extra_sel.empty())
    {
        debug_log("applying extra selection: " + extra_sel);
        e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(extra_sel);
        e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(extra_sel);
        if (p_data != nullptr)
            p_data->selection.nominal.node = p_data->selection.nominal.node.Filter(extra_sel);
    }

    Plotter plotter;
    auto &opt = plotter.options();
    opt.use_log_y = use_logy;
    opt.legend_on_top = true;
    opt.annotate_numbers = true;
    opt.overlay_signal = true;
    opt.show_ratio = include_data;
    opt.show_ratio_band = include_data;
    opt.adaptive_binning = true;
    opt.adaptive_min_sumw = 150.0;
    opt.adaptive_max_relerr = 0.20;
    opt.adaptive_fold_overflow = true;
    opt.signal_channels = Channels::signal_keys();
    opt.y_title = "Events/bin";
    opt.run_numbers = {"1"};
    opt.image_format = "pdf";

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();

    debug_log("plot options configured: include_data=" + std::string(include_data ? "true" : "false") +
              ", use_logy=" + std::string(use_logy ? "true" : "false"));

    struct DynamicAxis
    {
        int nbins = 50;
        double xmin = 0.0;
        double xmax = 1.0;
    };

    const auto build_dynamic_axis = [&](const std::string &expr,
                                        int fallback_nbins,
                                        double fallback_xmin,
                                        double fallback_xmax) {
        DynamicAxis out;
        out.nbins = fallback_nbins;
        out.xmin = fallback_xmin;
        out.xmax = fallback_xmax;

        double global_min = std::numeric_limits<double>::infinity();
        double global_max = -std::numeric_limits<double>::infinity();

        const auto update_minmax = [&](ROOT::RDF::RNode node) {
            const auto n_evt = node.Count().GetValue();
            if (n_evt == 0)
            {
                return;
            }

            const double local_min = static_cast<double>(node.Min(expr).GetValue());
            const double local_max = static_cast<double>(node.Max(expr).GetValue());
            global_min = std::min(global_min, local_min);
            global_max = std::max(global_max, local_max);
        };

        update_minmax(e_mc.selection.nominal.node);
        update_minmax(e_ext.selection.nominal.node);
        if (p_data != nullptr)
        {
            update_minmax(p_data->selection.nominal.node);
        }

        if (!std::isfinite(global_min) || !std::isfinite(global_max) || global_max <= global_min)
        {
            return out;
        }

        const double fallback_span = std::max(1.0, fallback_xmax - fallback_xmin);
        const double nominal_fine_width = fallback_span / static_cast<double>(fallback_nbins);
        const double span = std::max(nominal_fine_width, global_max - global_min);

        int dynamic_nbins = static_cast<int>(std::ceil(span / nominal_fine_width)) + 2;
        dynamic_nbins = std::max(12, dynamic_nbins);
        out.nbins = dynamic_nbins;

        const double padded_width = span / static_cast<double>(std::max(1, dynamic_nbins - 2));
        out.xmin = global_min - padded_width;
        out.xmax = global_max + padded_width;

        return out;
    };

    const auto draw_one = [&](const std::string &expr,
                              int nbins,
                              double xmin,
                              double xmax,
                              const std::string &x_title,
                              bool add_leading_empty_bin = false) {
        DynamicAxis axis = build_dynamic_axis(expr, nbins, xmin, xmax);

        std::string draw_expr = expr;
        if (add_leading_empty_bin && axis.nbins > 0)
        {
            const double bin_width = (axis.xmax - axis.xmin) / static_cast<double>(axis.nbins);
            axis.xmax += 2.0 * bin_width;
            axis.nbins += 2;
            draw_expr = "(" + expr + ") + " + std::to_string(bin_width);
        }

        opt.x_title = x_title.empty() ? expr : x_title;
        debug_log("drawing start: expr=" + expr +
                  ", draw_expr=" + draw_expr +
                  ", x_title=" + opt.x_title +
                  ", nbins=" + std::to_string(axis.nbins) +
                  ", xmin=" + std::to_string(axis.xmin) +
                  ", xmax=" + std::to_string(axis.xmax));

        const std::string weight = mc_weight.empty() ? "w_nominal" : mc_weight;
        TH1DModel spec = make_spec(expr, axis.nbins, axis.xmin, axis.xmax, weight);
        spec.expr = draw_expr;
        spec.sel = Preset::Empty;

        if (include_data)
        {
            plotter.draw_stack(spec, mc, data);
        }
        else
        {
            plotter.draw_stack(spec, mc);
        }

        debug_log("drawing done: expr=" + expr);
    };

    const int nbins = 50;
    draw_one("reco_neutrino_vertex_sce_z", nbins, -50.0, 1100.0, "Reco SCE neutrino vertex z [cm]", true);
    draw_one("reco_neutrino_vertex_sce_x", nbins, -50.0, 300.0, "Reco SCE neutrino vertex x [cm]", true);
    draw_one("reco_neutrino_vertex_sce_y", nbins, -180.0, 180.0, "Reco SCE neutrino vertex y [cm]", true);

    debug_log("completed all draw calls");
    return 0;
}


int plotStackedHistRecoSceNeutrinoVertex(const std::string &samples_tsv = "",
                              const std::string &extra_sel = "true",
                              bool use_logy = false,
                              bool include_data = false)
{
    const std::string mc_weight = "w_nominal";

    return plot_stacked_hist_impl(samples_tsv,
                                  mc_weight,
                                  extra_sel,
                                  use_logy,
                                  include_data);
}
