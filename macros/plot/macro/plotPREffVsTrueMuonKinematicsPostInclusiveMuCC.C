// plot/macro/plotPREffVsTrueMuonKinematicsPostInclusiveMuCC.C
//
// Plot pattern-recognition "assignment efficiency" vs true muon kinematics after inclusive #nu_{#mu} CC selection.
//
// Pattern-recognition efficiency definition (default), matching the paper:
//   - Additional event-level preselection:
//       * pass the inclusive #nu_{#mu} CC selection (sel_triggered_muon)
//
//   - Denominator (eligible truth events):
//       * nu_mu / anti-nu_mu CC   (is_nu_mu_cc)
//       * true neutrino vertex in fiducial volume (nu_vtx_in_fv)
//       * an in-event truth Lambda candidate (lam_pdg==3122) with p/pi daughters
//         (here approximated by requiring finite positive p_p and pi_p)
//
//   - Numerator (pass):
//       * a one-to-one assignment exists (pr_valid_assignment)
//       * for each (mu, p, pi): completeness > 0.1 and purity > 0.5
//
// Run with:
//   ./heron macro plotPREffVsTrueMuonKinematicsPostInclusiveMuCC.C
//   ./heron macro plotPREffVsTrueMuonKinematicsPostInclusiveMuCC.C 'plotPREffVsTrueMuonKinematicsPostInclusiveMuCC("./scratch/out/event_list_myana.root")'
//   ./heron macro plotPREffVsTrueMuonKinematicsPostInclusiveMuCC.C 'plotPREffVsTrueMuonKinematicsPostInclusiveMuCC("./scratch/out/event_list_myana.root","sel_triggered_muon && sel_reco_fv")'
//   ./heron macro plotPREffVsTrueMuonKinematicsPostInclusiveMuCC.C 'plotPREffVsTrueMuonKinematicsPostInclusiveMuCC("./scratch/out/event_list_myana.root","sel_triggered_muon","pr_valid_assignment && pr_mu_completeness>0.2 && pr_mu_purity>0.7 && pr_p_completeness>0.2 && pr_p_purity>0.7 && pr_pi_completeness>0.2 && pr_pi_purity>0.7")'
//
// Output:
//   Saves one plot per x-variable to:
//     $HERON_PLOT_DIR (default: ./scratch/plots)
//   with format:
//     $HERON_PLOT_FORMAT (default: pdf)

#include <algorithm>
#include <cstdlib>
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
#include "AnalysisContext.hh"
#include "ExecutionPolicy.hh"

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

struct VarSpec
{
    std::string expr;
    int nbins = 30;
    double xmin = 0.0;
    double xmax = 1.0;
    std::string x_title;
};

bool env_truthy(const char *v)
{
    if (v == nullptr)
        return false;
    const std::string s(v);
    return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
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

    std::cerr << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] missing required columns for " << label << ":\n";
    for (const auto &name : missing)
        std::cerr << "  - " << name << "\n";
    return 1;
}

} // namespace

int plotPREffVsTrueMuonKinematicsPostInclusiveMuCC(const std::string &samples_tsv = "",
                                                   const std::string &extra_sel = "sel_triggered_muon",
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
    const bool use_imt = env_truthy(std::getenv("HERON_ENABLE_IMT"));
    const ExecutionPolicy policy{.enableImplicitMT = use_imt};
    heron::AnalysisContext<ExecutionPolicy> context(policy, __func__);
    context.apply_runtime(__func__);

    std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] implicit MT "
              << (context.implicit_mt_enabled() ? "enabled via HERON_ENABLE_IMT" : "disabled (set HERON_ENABLE_IMT=1 to enable)")
              << "\n";

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] opening event list\n";
    EventListIO el(list_path);
    std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] building RDF from tree='" << el.event_tree() << "'\n";
    ROOT::RDataFrame rdf = el.rdf();
    std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] loaded sample refs: " << el.sample_refs().size() << "\n";

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    if (mask_mc)
        std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] MC-like mask size=" << mask_mc->size() << "\n";
    if (mask_ext)
        std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] EXT mask size=" << mask_ext->size() << "\n";

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
    std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] available columns=" << columns.size() << "\n";

    const std::vector<std::string> required_for_mask = {
        "sample_id"
    };
    if (require_columns(columns, required_for_mask, "sample filtering") != 0)
        return 1;

    const std::vector<std::string> required_for_denom_and_pass = {
        "sel_triggered_muon",
        "is_nu_mu_cc",
        "nu_vtx_in_fv",
        "lam_pdg",
        "p_p",
        "pi_p",
        "pr_valid_assignment",
        "pr_mu_completeness",
        "pr_mu_purity",
        "pr_p_completeness",
        "pr_p_purity",
        "pr_pi_completeness",
        "pr_pi_purity"
    };
    if (require_columns(columns, required_for_denom_and_pass, "denominator/pass selections") != 0)
        return 1;

    const std::vector<VarSpec> vars = {
        {"mu_p",      30, 0.0, 5.0, "#mu |p| [GeV/c]"},
        {"mu_theta",  30, 0.0, 3.2, "#theta(#nu,#mu) [rad]"},
    };

    int rc = 0;
    int plotted = 0;
    for (const auto &v : vars)
    {
        std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] processing variable='" << v.expr << "'\n";
        if (columns.find(v.expr) == columns.end())
        {
            std::cerr << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] skipping variable '" << v.expr
                      << "' (column not present in input tree)\n";
            rc = std::max(rc, 2);
            continue;
        }

        TH1DModel spec = make_spec(v.expr, v.nbins, v.xmin, v.xmax, "");

        EfficiencyPlot::Config cfg;
        cfg.x_title = v.x_title;
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
        opt.out_dir = context.outputs().plot_dir;
        opt.image_format.clear();

        EfficiencyPlot eff(spec, opt, cfg);
        std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] computing efficiency for '" << v.expr << "'\n";
        const int one = eff.compute(node_mc, denom_sel, pass_sel, extra_sel);
        std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] compute rc=" << one
                  << " denom=" << eff.denom_entries()
                  << " pass=" << eff.pass_entries() << "\n";
        if (one == 0)
        {
            const std::string tag = Plotter::sanitise(v.expr);
            const int drc = eff.draw_and_save("pr_eff_post_inclusive_mu_cc_true_mu_vs_" + tag);
            std::cout << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] draw rc=" << drc
                      << " variable='" << v.expr << "'\n";
            rc = std::max(rc, drc);
            ++plotted;
        }
        else
        {
            rc = std::max(rc, one);
        }
    }

    if (plotted == 0)
    {
        std::cerr << "[plotPREffVsTrueMuonKinematicsPostInclusiveMuCC] no plots were produced; check input column availability\n";
        rc = std::max(rc, 1);
    }

    return rc;
}
