// plot/macro/plotStackedHistRecoSceNeutrinoVertex.C
//
// Run with:
//   ./heron macro plotStackedHistRecoSceNeutrinoVertex.C
//   ./heron macro plotStackedHistRecoSceNeutrinoVertex.C 'plotStackedHistRecoSceNeutrinoVertex("/path/to/samples.tsv","true",false)'
//
// Notes:
//   - This macro loads aggregated samples (samples.tsv -> SampleIO -> original analysis tree)
//   - It runs your analysis column derivations so that "analysis_channels" exists for stacking.
//   - The stack is grouped by "analysis_channels"; expr controls the x-axis variable only.
//   - MC yields are scaled by w_nominal unless an alternative weight is provided.
//   - Output dir/format follow PlotEnv defaults (HERON_PLOT_DIR / HERON_PLOT_FORMAT).
//   - Default input uses the generated event list (event_list_<analysis>.root).

#include <cstdlib>
#include <iostream>
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
#include "MacroGuard.hh"


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
    const char *env = std::getenv("HERON_DEBUG_PLOT_STACK");
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

    const std::string software_trigger_gate_sel = "software_trigger > 0";
    const std::string reco_neutrino_slice_sel = "sel_slice";
    const std::string gate_sel = "(" + software_trigger_gate_sel + ") && (" + reco_neutrino_slice_sel + ")";
    const std::string combined_sel = extra_sel.empty()
                                         ? gate_sel
                                         : "(" + gate_sel + ") && (" + extra_sel + ")";
    debug_log("applying selection (software trigger + neutrino slice gates): " + combined_sel);
    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(combined_sel);
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(combined_sel);
    if (p_data != nullptr)
        p_data->selection.nominal.node = p_data->selection.nominal.node.Filter(combined_sel);

    Plotter plotter;
    auto &opt = plotter.options();
    opt.use_log_y = use_logy;
    opt.legend_on_top = true;
    opt.annotate_numbers = true;
    opt.overlay_signal = true;
    opt.show_ratio = include_data;
    opt.show_ratio_band = include_data;
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

    constexpr double reco_fv_min_x = 5.0;
    constexpr double reco_fv_max_x = 251.0;
    constexpr double reco_fv_min_y = -110.0;
    constexpr double reco_fv_max_y = 110.0;
    constexpr double reco_fv_min_z = 20.0;
    constexpr double reco_fv_max_z = 986.0;
    constexpr double reco_gap_min_z = 675.0;
    constexpr double reco_gap_max_z = 775.0;

    const auto draw_one = [&](const std::string &expr,
                              int nbins,
                              double xmin,
                              double xmax,
                              const std::string &x_title,
                              bool add_leading_empty_bin = false) {
        int axis_nbins = nbins;
        double axis_xmin = xmin;
        double axis_xmax = xmax;

        std::string draw_expr = expr;
        if (add_leading_empty_bin && axis_nbins > 0)
        {
            const double bin_width = (axis_xmax - axis_xmin) / static_cast<double>(axis_nbins);
            axis_xmax += 2.0 * bin_width;
            axis_nbins += 2;
            draw_expr = "(" + expr + ") + " + std::to_string(bin_width);
        }

        opt.x_title = x_title.empty() ? expr : x_title;

        opt.show_cuts = true;
        if (expr == "reco_neutrino_vertex_sce_x")
        {
            opt.cuts = {{reco_fv_min_x, CutDir::GreaterThan},
                        {reco_fv_max_x, CutDir::LessThan}};
        }
        else if (expr == "reco_neutrino_vertex_sce_y")
        {
            opt.cuts = {{reco_fv_min_y, CutDir::GreaterThan},
                        {reco_fv_max_y, CutDir::LessThan}};
        }
        else if (expr == "reco_neutrino_vertex_sce_z")
        {
            opt.cuts = {{reco_fv_min_z, CutDir::GreaterThan},
                        {reco_gap_min_z, CutDir::LessThan},
                        {reco_gap_max_z, CutDir::GreaterThan},
                        {reco_fv_max_z, CutDir::LessThan}};
        }
        else
        {
            opt.show_cuts = false;
            opt.cuts.clear();
        }

        debug_log("drawing start: expr=" + expr +
                  ", draw_expr=" + draw_expr +
                  ", x_title=" + opt.x_title +
                  ", axis_nbins=" + std::to_string(axis_nbins) +
                  ", axis_xmin=" + std::to_string(axis_xmin) +
                  ", axis_xmax=" + std::to_string(axis_xmax));

        const std::string weight = mc_weight.empty() ? "w_nominal" : mc_weight;
        TH1DModel spec = make_spec(expr, axis_nbins, axis_xmin, axis_xmax, weight);
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

    draw_one("reco_neutrino_vertex_sce_z", 230, -50.0, 1100.0, "Reco SCE neutrino vertex z [cm]", true);
    draw_one("reco_neutrino_vertex_sce_x", 70, -50.0, 300.0, "Reco SCE neutrino vertex x [cm]", true);
    draw_one("reco_neutrino_vertex_sce_y", 72, -180.0, 180.0, "Reco SCE neutrino vertex y [cm]", true);

    debug_log("completed all draw calls");
    return 0;
}


int plotStackedHistRecoSceNeutrinoVertex(const std::string &samples_tsv = "",
                              const std::string &extra_sel = "true",
                              bool use_logy = false,
                              bool include_data = false)
{
  return heron::macro::run_with_guard("plotStackedHistRecoSceNeutrinoVertex", [&]() -> int {
    const std::string mc_weight = "w_nominal";

    return plot_stacked_hist_impl(samples_tsv,
                                  mc_weight,
                                  extra_sel,
                                  use_logy,
                                  include_data);

  });
}
