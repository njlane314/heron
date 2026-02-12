// plot/macro/plotParticleLevelMipnessSelectionStacked.C
//
// Particle-level stacked histogram of track MIPness, using the cumulative
// selection up to (and including) the reco-fiducial stage by default.
//
// Run with:
//   ./nuxsec macro plotParticleLevelMipnessSelectionStacked.C
//   ./nuxsec macro plotParticleLevelMipnessSelectionStacked.C \
//     'plotParticleLevelMipnessSelectionStacked("./scratch/out/event_list_myana.root")'
//   ./nuxsec macro plotParticleLevelMipnessSelectionStacked.C \
//     'plotParticleLevelMipnessSelectionStacked("", "w_nominal", "", false, true, "sel_trigger && sel_triggered_slice && sel_reco_fv && sel_triggered_muon")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
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


} // namespace

int plotParticleLevelMipnessSelectionStacked(const std::string &event_list_path = "",
                                             const std::string &mc_weight = "w_nominal",
                                             const std::string &extra_sel = "",
                                             bool use_logy = false,
                                             bool include_data = true,
                                             const std::string &selection_upto = "sel_trigger && sel_triggered_slice && sel_reco_fv",
                                             const std::string &particle_pdg_branch = "backtracked_pdg_codes",
                                             double squash_midpoint = 1.0,
                                             double squash_softness = 0.2)
{
    ROOT::EnableImplicitMT();

    const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;

    if (!looks_like_event_list_root(input_path))
    {
        std::cerr << "[plotParticleLevelMipnessSelectionStacked] input is not an event-list root file: " << input_path << "\n";
        return 1;
    }

    EventListIO el(input_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

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

    ProcessorEntry rec_mc;
    rec_mc.source = Type::kMC;

    ProcessorEntry rec_ext;
    rec_ext.source = Type::kExt;

    ProcessorEntry rec_data;
    rec_data.source = Type::kData;

    auto has_mipness_column = [](ROOT::RDF::RNode node) {
        const auto names = node.GetColumnNames();
        return std::find(names.begin(), names.end(), "track_mipness_median_plane_score") != names.end();
    };

    if (!has_mipness_column(node_mc) || !has_mipness_column(node_ext) || (include_data && !has_mipness_column(node_data)))
    {
        std::cerr << "[plotParticleLevelMipnessSelectionStacked] missing required column: "
                  << "track_mipness_median_plane_score\n";
        return 1;
    }

    if (squash_softness <= 0.0)
    {
        std::cerr << "[plotParticleLevelMipnessSelectionStacked] squash_softness must be > 0, got "
                  << squash_softness << "\n";
        return 1;
    }

    std::vector<Entry> entries;
    entries.reserve(include_data ? 3 : 2);

    std::vector<const Entry *> mc;
    std::vector<const Entry *> data;

    entries.emplace_back(make_entry(std::move(node_mc), rec_mc));
    Entry &e_mc = entries.back();
    mc.push_back(&e_mc);

    entries.emplace_back(make_entry(std::move(node_ext), rec_ext));
    Entry &e_ext = entries.back();
    mc.push_back(&e_ext);

    Entry *p_data = nullptr;
    if (include_data)
    {
        entries.emplace_back(make_entry(std::move(node_data), rec_data));
        p_data = &entries.back();
        data.push_back(p_data);
    }

    const std::string stage_sel = selection_upto.empty() ? std::string("true") : selection_upto;
    const std::string combined_sel = extra_sel.empty()
                                         ? "(" + stage_sel + ")"
                                         : "(" + stage_sel + ") && (" + extra_sel + ")";

    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(combined_sel);
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(combined_sel);
    if (p_data != nullptr)
        p_data->selection.nominal.node = p_data->selection.nominal.node.Filter(combined_sel);

    std::ostringstream squashed_expr;
    squashed_expr << "1.0/(1.0 + exp(("
                  << squash_midpoint
                  << " - track_mipness_median_plane_score)/"
                  << squash_softness
                  << "))";
    const std::string squashed_col = "track_mipness_squashed_score";

    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Define(squashed_col, squashed_expr.str());
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Define(squashed_col, squashed_expr.str());
    if (p_data != nullptr)
        p_data->selection.nominal.node = p_data->selection.nominal.node.Define(squashed_col, squashed_expr.str());

    Plotter plotter;
    auto &opt = plotter.options();
    opt.use_log_y = use_logy;
    opt.legend_on_top = true;
    opt.annotate_numbers = true;
    opt.overlay_signal = false;
    opt.show_ratio = include_data;
    opt.show_ratio_band = include_data;
    opt.signal_channels = Channels::signal_keys();
    opt.x_title = "Track MIPness score";
    opt.y_title = "Particles";
    opt.analysis_region_label = "Selection up to MIPness stage";
    opt.image_format = "pdf";
    opt.particle_level = true;
    opt.particle_pdg_branch = particle_pdg_branch;

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();
    opt.run_numbers = {"1"};

    TH1DModel spec;
    spec.id = "particle_level_track_mipness_selection_upto";
    spec.name = "particle_level_track_mipness_selection_upto";
    spec.title = "Track MIPness (particle-level stack)";
    spec.expr = "track_mipness_median_plane_score";
    spec.weight = mc_weight;
    spec.nbins = 50;
    spec.xmin = 0.0;
    spec.xmax = 2.0;
    spec.sel = Preset::Empty;

    std::cout << "[plotParticleLevelMipnessSelectionStacked] selection=" << combined_sel
              << ", include_data=" << (include_data ? "true" : "false")
              << ", use_logy=" << (use_logy ? "true" : "false")
              << ", squash_midpoint=" << squash_midpoint
              << ", squash_softness=" << squash_softness
              << ", particle_pdg_branch=" << opt.particle_pdg_branch << "\n";

    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);

    std::cout << "[plotParticleLevelMipnessSelectionStacked] wrote "
              << plotter.options().out_dir << "/" << spec.id << "."
              << plotter.options().image_format << "\n";

    TH1DModel spec_squashed = spec;
    spec_squashed.id = "particle_level_track_mipness_selection_upto_squashed";
    spec_squashed.name = "particle_level_track_mipness_selection_upto_squashed";
    spec_squashed.title = "Squashed track MIPness (particle-level stack)";
    spec_squashed.expr = squashed_col;
    spec_squashed.nbins = 50;
    spec_squashed.xmin = 0.0;
    spec_squashed.xmax = 1.0;
    opt.x_title = "Squashed track MIPness score";

    if (include_data)
        plotter.draw_stack(spec_squashed, mc, data);
    else
        plotter.draw_stack(spec_squashed, mc);

    std::cout << "[plotParticleLevelMipnessSelectionStacked] wrote "
              << plotter.options().out_dir << "/" << spec_squashed.id << "."
              << plotter.options().image_format << "\n";

    return 0;
}
