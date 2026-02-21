// plot/macro/plotSoftwareTriggerGatePassFail.C
//
// Stacked pass/fail plot for the software-trigger gate from EventListIO output.
//
// The plotted variable is:
//   software_trigger_gate_pass = (software_trigger > 0 ? 1 : 0)
//
// Run with:
//   ./heron macro plotSoftwareTriggerGatePassFail.C
//   ./heron macro plotSoftwareTriggerGatePassFail.C \
//     'plotSoftwareTriggerGatePassFail("./scratch/out/event_list_myana.root")'
//   ./heron macro plotSoftwareTriggerGatePassFail.C \
//     'plotSoftwareTriggerGatePassFail("./scratch/out/event_list_myana.root","sel_slice",true,true)'

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include <TFile.h>

#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
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

std::string getenv_or(const char *key, const std::string &fallback)
{
    const char *v = std::getenv(key);
    if (!v || std::string(v).empty())
        return fallback;
    return std::string(v);
}

ROOT::RDF::RNode filter_by_mask(ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask)
{
    if (!mask)
        return n;

    return n.Filter(
        [mask](int sid) {
            return sid >= 0
                   && sid < static_cast<int>(mask->size())
                   && (*mask)[static_cast<size_t>(sid)];
        },
        {"sample_id"});
}

} // namespace

int plotSoftwareTriggerGatePassFail(const std::string &input = "",
                                    const std::string &extra_sel = "true",
                                    bool use_logy = false,
                                    bool include_data = true)
{
  return heron::macro::run_with_guard("plotSoftwareTriggerGatePassFail", [&]() -> int {
    const std::string list_path = input.empty() ? default_event_list_root() : input;
    std::cout << "[plotSoftwareTriggerGatePassFail] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotSoftwareTriggerGatePassFail] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

    ROOT::RDF::RNode node_ext = filter_by_mask(rdf, mask_ext);
    ROOT::RDF::RNode node_mc = filter_by_mask(rdf, mask_mc);
    if (mask_ext)
    {
        node_mc = node_mc.Filter(
            [mask_ext](int sid) {
                return !(sid >= 0
                         && sid < static_cast<int>(mask_ext->size())
                         && (*mask_ext)[static_cast<size_t>(sid)]);
            },
            {"sample_id"});
    }
    ROOT::RDF::RNode node_data = filter_by_mask(rdf, mask_data);

    std::vector<Entry> entries;
    entries.reserve(include_data ? 3 : 2);

    std::vector<const Entry *> mc;
    std::vector<const Entry *> data;

    ProcessorEntry rec_mc;
    rec_mc.source = Type::kMC;
    entries.emplace_back(make_entry(std::move(node_mc), rec_mc));
    Entry &e_mc = entries.back();
    mc.push_back(&e_mc);

    ProcessorEntry rec_ext;
    rec_ext.source = Type::kExt;
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

    const std::string selection = (extra_sel.empty() ? "true" : extra_sel);
    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(selection);
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(selection);
    if (p_data)
        p_data->selection.nominal.node = p_data->selection.nominal.node.Filter(selection);

    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Define(
        "software_trigger_gate_pass", [](int sw) { return sw > 0 ? 1.0 : 0.0; }, {"software_trigger"});
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Define(
        "software_trigger_gate_pass", [](int sw) { return sw > 0 ? 1.0 : 0.0; }, {"software_trigger"});
    if (p_data)
    {
        p_data->selection.nominal.node = p_data->selection.nominal.node.Define(
            "software_trigger_gate_pass", [](int sw) { return sw > 0 ? 1.0 : 0.0; }, {"software_trigger"});
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
    opt.show_cuts = true;
    opt.cuts = {{0.5, CutDir::GreaterThan}};
    opt.x_title = "Slice Selected";
    opt.y_title = "Events";
    opt.analysis_region_label = "Software Trigger Gate";
    opt.image_format = getenv_or("HERON_PLOT_FORMAT", "pdf");

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();
    opt.run_numbers = {"1"};

    TH1DModel spec;
    spec.id = "software_trigger_gate_pass_fail";
    spec.name = "software_trigger_gate_pass_fail";
    spec.title = "Software Trigger Gate";
    spec.expr = "software_trigger_gate_pass";
    spec.weight = "w_nominal";
    spec.nbins = 2;
    spec.xmin = -0.5;
    spec.xmax = 1.5;
    spec.sel = Preset::Empty;

    std::cout << "[plotSoftwareTriggerGatePassFail] selection=" << selection
              << ", include_data=" << (include_data ? "true" : "false")
              << ", use_logy=" << (use_logy ? "true" : "false") << "\n";

    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);

    std::cout << "[plotSoftwareTriggerGatePassFail] wrote "
              << plotter.options().out_dir << "/" << spec.id << "."
              << plotter.options().image_format << "\n";

    return 0;

  });
}
