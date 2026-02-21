// plot/macro/plotPREffLambdaSignalVsDecayDistance.C
//
// Plot pattern-recognition assignment efficiency for true #Lambda signal
// versus true decay distance from the neutrino vertex, after neutrino-slice selection.
//
// Default event-level preselection:
//   sel_slice
//
// Pattern-recognition efficiency definition (default):
//   - Denominator (eligible truth events):
//       is_nu_mu_cc && nu_vtx_in_fv && (lam_pdg==3122) && (p_p>0.0) && (pi_p>0.0)
//   - Numerator (pass):
//       pr_valid_assignment
//       && pr_mu_completeness>0.1 && pr_mu_purity>0.5
//       && pr_p_completeness>0.1 && pr_p_purity>0.5
//       && pr_pi_completeness>0.1 && pr_pi_purity>0.5
//
// Run with:
//   ./heron macro plotPREffLambdaSignalVsDecayDistance.C
//   ./heron macro plotPREffLambdaSignalVsDecayDistance.C \
//     'plotPREffLambdaSignalVsDecayDistance("./scratch/out/event_list_myana.root")'

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TEfficiency.h>
#include <TFile.h>

#include "EfficiencyPlot.hh"
#include "EventListIO.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "include/MacroGuard.hh"

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

int require_columns(const std::unordered_set<std::string> &columns,
                    const std::vector<std::string> &required,
                    const std::string &label)
{
    std::vector<std::string> missing;
    missing.reserve(required.size());

    for (const auto &name : required)
    {
        if (columns.find(name) == columns.end())
            missing.push_back(name);
    }

    if (missing.empty())
        return 0;

    std::cerr << "[plotPREffLambdaSignalVsDecayDistance] missing required columns for " << label << ":\n";
    for (const auto &name : missing)
        std::cerr << "  - " << name << "\n";

    return 1;
}

} // namespace

int plotPREffLambdaSignalVsDecayDistance(
    const std::string &samples_tsv = "",
    const std::string &extra_sel = "sel_slice",
    const std::string &pass_sel =
        "pr_valid_assignment"
        " && pr_mu_completeness>0.1 && pr_mu_purity>0.5"
        " && pr_p_completeness>0.1 && pr_p_purity>0.5"
        " && pr_pi_completeness>0.1 && pr_pi_purity>0.5",
    const std::string &denom_sel =
        "is_nu_mu_cc"
        " && nu_vtx_in_fv"
        " && (lam_pdg==3122)"
        " && (p_p>0.0)"
        " && (pi_p>0.0)")
{
  return heron::macro::run_with_guard("plotPREffLambdaSignalVsDecayDistance", [&]() -> int {
    ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotPREffLambdaSignalVsDecayDistance] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotPREffLambdaSignalVsDecayDistance] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();

    auto filter_by_mask = [](ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask) -> ROOT::RDF::RNode {
        if (!mask)
            return n;
        return n.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<size_t>(sid)];
            },
            {"sample_id"});
    };

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

    const auto col_names = node_mc.GetColumnNames();
    const std::unordered_set<std::string> columns(col_names.begin(), col_names.end());

    const std::vector<std::string> required_cols = {
        "sample_id",
        "sel_slice",
        "is_nu_mu_cc",
        "nu_vtx_in_fv",
        "lam_pdg",
        "lam_decay_sep",
        "p_p",
        "pi_p",
        "pr_valid_assignment",
        "pr_mu_completeness",
        "pr_mu_purity",
        "pr_p_completeness",
        "pr_p_purity",
        "pr_pi_completeness",
        "pr_pi_purity"};
    if (require_columns(columns, required_cols, "selections and plotting") != 0)
        return 1;

    TH1DModel spec = make_spec("lam_decay_sep", 30, 0.0, 20.0, "");

    EfficiencyPlot::Config cfg;
    cfg.x_title = "#Lambda decay sep. to #nu vtx [cm]";
    cfg.y_counts_title = "Events";
    cfg.y_eff_title = "Efficiency";
    cfg.legend_total = "Total";
    cfg.legend_passed = "Passed";
    cfg.legend_eff = "Efficiency";
    cfg.stat = TEfficiency::kFCP;
    cfg.conf_level = 0.68;
    cfg.auto_x_range = false;
    cfg.draw_distributions = true;

    Options opt;
    opt.out_dir.clear();
    opt.image_format.clear();

    EfficiencyPlot eff(spec, opt, cfg);
    const int rc = eff.compute(node_mc, denom_sel, pass_sel, extra_sel);
    if (rc != 0)
        return rc;

    return eff.draw_and_save("pr_eff_lambda_signal_vs_decay_distance_post_neutrino_slice");

  });
}
