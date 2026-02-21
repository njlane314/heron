// plot/macro/plotStackedHistRecoSceNeutrinoVertexEmptySelection.C
//
// Run with:
//   ./heron macro plotStackedHistRecoSceNeutrinoVertexEmptySelection.C
//   ./heron macro plotStackedHistRecoSceNeutrinoVertexEmptySelection.C 'plotStackedHistRecoSceNeutrinoVertexEmptySelection("/path/to/event_list.root",false,false)'
//
// Notes:
//   - This macro loads the event list ROOT file and stacks channels using analysis_channels.
//   - The plots use the empty histogram preset after software-trigger and slice gates.
//   - The x-axis ranges match the standard true-vertex projection plots.
//   - Binning is fixed and uniform in each projection.

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>

#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "MacroGuard.hh"

using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &path)
{
    const auto n = path.size();
    if (n < 5 || path.substr(n - 5) != ".root")
    {
        return false;
    }

    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie())
    {
        return false;
    }

    const bool has_refs = (file->Get("sample_refs") != nullptr);
    const bool has_events_tree = (file->Get("events") != nullptr);
    const bool has_event_tree_key = (file->Get("event_tree") != nullptr);

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

    std::cout << "[plotStackedHistRecoSceNeutrinoVertexEmptySelection][debug] " << msg << "\n";
    std::cout.flush();
}

bool implicit_mt_enabled()
{
    const char *env = std::getenv("HERON_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}
} // namespace

int plotStackedHistRecoSceNeutrinoVertexEmptySelection(const std::string &event_list_root = "",
                                                       bool use_logy = false,
                                                       bool include_data = false)
{
  return heron::macro::run_with_guard("plotStackedHistRecoSceNeutrinoVertexEmptySelection", [&]() -> int {
    if (implicit_mt_enabled())
    {
        ROOT::EnableImplicitMT();
        debug_log("ROOT implicit MT enabled (HERON_PLOT_IMT != 0)");
    }
    else
    {
        debug_log("ROOT implicit MT disabled (set HERON_PLOT_IMT=1 to enable)");
    }

    const std::string list_path = event_list_root.empty() ? default_event_list_root() : event_list_root;
    std::cout << "[plotStackedHistRecoSceNeutrinoVertexEmptySelection] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotStackedHistRecoSceNeutrinoVertexEmptySelection] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

    auto filter_by_mask = [](ROOT::RDF::RNode node, std::shared_ptr<const std::vector<char>> mask) {
        return node.Filter(
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
    Entry *e_data = nullptr;

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

    if (include_data)
    {
        ProcessorEntry rec_data;
        rec_data.source = Type::kData;
        entries.emplace_back(make_entry(std::move(node_data), rec_data));
        e_data = &entries.back();
        data.push_back(e_data);
    }

    const std::string software_trigger_gate_sel = "software_trigger > 0";
    const std::string reco_neutrino_slice_sel = "sel_slice";
    const std::string combined_gate_sel = "(" + software_trigger_gate_sel + ") && (" + reco_neutrino_slice_sel + ")";
    debug_log("applying selection (software trigger + neutrino slice gates): " + combined_gate_sel);
    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(combined_gate_sel);
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(combined_gate_sel);
    if (include_data)
    {
        e_data->selection.nominal.node = e_data->selection.nominal.node.Filter(combined_gate_sel);
    }

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

    const auto draw_one = [&](const std::string &expr,
                              int nbins,
                              double xmin,
                              double xmax,
                              const std::string &x_title) {
        opt.x_title = x_title;

        TH1DModel spec = make_spec(expr, nbins, xmin, xmax, "w_nominal");
        spec.expr = expr;
        spec.sel = Preset::Empty;

        if (include_data)
        {
            plotter.draw_stack(spec, mc, data);
        }
        else
        {
            plotter.draw_stack(spec, mc);
        }
    };

    const int nbins = 50;
    draw_one("reco_neutrino_vertex_sce_x", nbins, -50.0, 300.0, "Reco SCE neutrino vertex x [cm]");
    draw_one("reco_neutrino_vertex_sce_y", nbins, -180.0, 180.0, "Reco SCE neutrino vertex y [cm]");
    draw_one("reco_neutrino_vertex_sce_z", nbins, -50.0, 1100.0, "Reco SCE neutrino vertex z [cm]");

    return 0;

  });
}
