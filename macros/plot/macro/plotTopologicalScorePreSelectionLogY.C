// plot/macro/plotTopologicalScorePreSelectionLogY.C
//
// Stacked topological-score distribution before any selection (log-y by default).
//
// Run with:
//   ./heron macro plotTopologicalScorePreSelectionLogY.C
//   ./heron macro plotTopologicalScorePreSelectionLogY.C \
//     'plotTopologicalScorePreSelectionLogY("./scratch/out/event_list_myana.root")'
//   ./heron macro plotTopologicalScorePreSelectionLogY.C \
//     'plotTopologicalScorePreSelectionLogY("./scratch/out/event_list_myana.root","true",true,true)'

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
#include "SampleCLI.hh"
#include "AnalysisContext.hh"
#include "ExecutionPolicy.hh"

using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &path)
{
    const auto n = path.size();
    if (n < 5 || path.substr(n - 5) != ".root")
        return false;

    std::unique_ptr<TFile> input_file(TFile::Open(path.c_str(), "READ"));
    if (!input_file || input_file->IsZombie())
        return false;

    const bool has_refs = (input_file->Get("sample_refs") != nullptr);
    const bool has_events_tree = (input_file->Get("events") != nullptr);
    const bool has_event_tree_key = (input_file->Get("event_tree") != nullptr);
    return has_refs && (has_events_tree || has_event_tree_key);
}

}

int plotTopologicalScorePreSelectionLogY(const std::string &samples_tsv = "",
                                         const std::string &extra_sel = "true",
                                         bool use_logy = true,
                                         bool include_data = false)
{
    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotTopologicalScorePreSelectionLogY] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotTopologicalScorePreSelectionLogY] input is not an event list ROOT file: "
                  << list_path << "\n";
        return 2;
    }

    const ExecutionPolicy policy = ExecutionPolicy::from_env();
    heron::AnalysisContext<ExecutionPolicy> context(policy, __func__);
    context.apply_runtime(__func__);

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
                       && (*mask)[static_cast<std::size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base = rdf;
    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext);
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0
                                                && sid < static_cast<int>(mask_ext->size())
                                                && (*mask_ext)[static_cast<std::size_t>(sid)]);
                                   },
                                   {"sample_id"});
    ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

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

    if (!extra_sel.empty())
    {
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
    opt.signal_channels = Channels::signal_keys();
    opt.y_title = "Events";
    opt.run_numbers = {"1"};
    opt.image_format = "pdf";

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();

    TH1DModel spec = make_spec("topological_score", 80, 0.0, 1.0, "w_nominal");
    spec.sel = Preset::Empty;

    opt.x_title = "Topological score";

    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);

    std::cout << "[plotTopologicalScorePreSelectionLogY] done\n";
    return 0;
}
