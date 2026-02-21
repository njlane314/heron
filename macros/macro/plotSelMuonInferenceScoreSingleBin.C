// plot/macro/plotSelMuonInferenceScoreSingleBin.C
//
// Build a single-bin stacked histogram of events that pass:
//   sel_muon && inf_scores[0] > threshold
//
// Run with:
//   ./heron macro plotSelMuonInferenceScoreSingleBin.C
//   ./heron macro plotSelMuonInferenceScoreSingleBin.C 'plotSelMuonInferenceScoreSingleBin("./scratch/out/event_list_myana.root", 6.9, false, false)'

#include <cmath>
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

std::string make_selection_expr(double threshold)
{
    return "sel_muon && (inf_scores.size() > 0) && (inf_scores[0] > " + std::to_string(threshold) + ")";
}
} // namespace

int plotSelMuonInferenceScoreSingleBin(const std::string &event_list_root = "",
                                       double first_inference_score_cut = 6.9,
                                       bool use_logy = false,
                                       bool include_data = false)
{
  return heron::macro::run_with_guard("plotSelMuonInferenceScoreSingleBin", [&]() -> int {
    const std::string list_path = event_list_root.empty() ? default_event_list_root() : event_list_root;
    if (!looks_like_event_list_root(list_path))
    {
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

    const std::string sel_expr = make_selection_expr(first_inference_score_cut);

    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(sel_expr)
                                                   .Define("single_bin_axis", []() { return 0.5; });
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(sel_expr)
                                                     .Define("single_bin_axis", []() { return 0.5; });
    if (include_data)
    {
        e_data->selection.nominal.node = e_data->selection.nominal.node.Filter(sel_expr)
                                                        .Define("single_bin_axis", []() { return 0.5; });
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
    opt.x_title = "sel_muon && inf_scores[0] > " + std::to_string(first_inference_score_cut);
    opt.run_numbers = {"1"};
    opt.image_format = "pdf";

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();

    TH1DModel spec = make_spec("sel_muon_infscore_single_bin", 1, 0.0, 1.0, "w_nominal");
    spec.expr = "single_bin_axis";
    spec.sel = Preset::Empty;

    if (include_data)
    {
        plotter.draw_stack(spec, mc, data);
    }
    else
    {
        plotter.draw_stack(spec, mc);
    }

    return 0;

  });
}
