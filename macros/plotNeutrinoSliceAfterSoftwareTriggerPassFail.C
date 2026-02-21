// plot/macro/plotNeutrinoSliceAfterSoftwareTriggerPassFail.C
//
// Stacked pass/fail plot for neutrino-slice selection after software-trigger gate.
//
// The plotted variable is:
//   neutrino_slice_pass = (sel_slice ? 1 : 0)
//
// Run with:
//   ./heron macro plotNeutrinoSliceAfterSoftwareTriggerPassFail.C
//   ./heron macro plotNeutrinoSliceAfterSoftwareTriggerPassFail.C \
//     'plotNeutrinoSliceAfterSoftwareTriggerPassFail("./scratch/out/event_list_myana.root")'
//   ./heron macro plotNeutrinoSliceAfterSoftwareTriggerPassFail.C \
//     'plotNeutrinoSliceAfterSoftwareTriggerPassFail("./scratch/out/event_list_myana.root","abs(reco_nu_vtx_x) < 250",true,true)'

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
#include "include/MacroGuard.hh"
#include "include/MacroEnv.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{


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

int plotNeutrinoSliceAfterSoftwareTriggerPassFail(const std::string &input = "",
                                                  const std::string &extra_sel = "true",
                                                  bool use_logy = false,
                                                  bool include_data = true)
{
  return heron::macro::run_with_guard("plotNeutrinoSliceAfterSoftwareTriggerPassFail", [&]() -> int {
    const std::string list_path = input.empty() ? default_event_list_root() : input;
    std::cout << "[plotNeutrinoSliceAfterSoftwareTriggerPassFail] input=" << list_path << "\n";

    if (!heron::macro::looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotNeutrinoSliceAfterSoftwareTriggerPassFail] input is not an event list root file: " << list_path << "\n";
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

    const std::string selection = "(" + (extra_sel.empty() ? std::string("true") : extra_sel) + ") && software_trigger > 0";
    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(selection);
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(selection);
    if (p_data)
        p_data->selection.nominal.node = p_data->selection.nominal.node.Filter(selection);

    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Define(
        "neutrino_slice_pass", [](bool sel) { return sel ? 1.0 : 0.0; }, {"sel_slice"});
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Define(
        "neutrino_slice_pass", [](bool sel) { return sel ? 1.0 : 0.0; }, {"sel_slice"});
    if (p_data)
    {
        p_data->selection.nominal.node = p_data->selection.nominal.node.Define(
            "neutrino_slice_pass", [](bool sel) { return sel ? 1.0 : 0.0; }, {"sel_slice"});
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
    opt.x_title = "Neutrino Slice Selected";
    opt.y_title = "Events";
    opt.analysis_region_label = "Software Trigger + Neutrino Slice";
    opt.image_format = heron::macro::getenv_or("HERON_PLOT_FORMAT", "pdf");

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();
    opt.run_numbers = {"1"};

    TH1DModel spec;
    spec.id = "neutrino_slice_after_software_trigger_pass_fail";
    spec.name = "neutrino_slice_after_software_trigger_pass_fail";
    spec.title = "Neutrino Slice after Software Trigger";
    spec.expr = "neutrino_slice_pass";
    spec.weight = "w_nominal";
    spec.nbins = 2;
    spec.xmin = -0.5;
    spec.xmax = 1.5;
    spec.sel = Preset::Empty;

    std::cout << "[plotNeutrinoSliceAfterSoftwareTriggerPassFail] selection=" << selection
              << ", include_data=" << (include_data ? "true" : "false")
              << ", use_logy=" << (use_logy ? "true" : "false") << "\n";

    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);

    std::cout << "[plotNeutrinoSliceAfterSoftwareTriggerPassFail] wrote "
              << plotter.options().out_dir << "/" << spec.id << "."
              << plotter.options().image_format << "\n";

    return 0;

  });
}
