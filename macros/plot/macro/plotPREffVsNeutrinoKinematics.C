// plot/macro/plotPREffVsNeutrinoKinematics.C
//
// Plot pattern-recognition "assignment efficiency" vs true neutrino kinematics.
//
// Pattern-recognition efficiency definition (default), matching the paper:
//   - Denominator (eligible truth events):
//       is_nu_mu_cc && nu_vtx_in_fv && (lam_pdg==3122) && (p_p>0) && (pi_p>0)
//   - Numerator (pass):
//       pr_valid_assignment AND, for each of (mu,p,pi),
//       completeness>0.1 AND purity>0.5
//
// Run with:
//   ./heron macro plotPREffVsNeutrinoKinematics.C
//   ./heron macro plotPREffVsNeutrinoKinematics.C 'plotPREffVsNeutrinoKinematics("./scratch/out/event_list_myana.root","sel_triggered_slice")'
//
// Output:
//   Saves one plot per x-variable to:
//     $HERON_PLOT_DIR (default: ./scratch/plots)
//   with format:
//     $HERON_PLOT_FORMAT (default: pdf)

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TEfficiency.h>

#include "EfficiencyPlot.hh"
#include "EventListIO.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"

using namespace nu;

namespace
{
struct VarSpec
{
    std::string expr;
    int nbins = 30;
    double xmin = 0.0;
    double xmax = 1.0;
    std::string x_title;
};

} // namespace

int plotPREffVsNeutrinoKinematics(const std::string &samples_tsv = "",
                                  const std::string &extra_sel = "true",
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
    ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotPREffVsNeutrinoKinematics] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotPREffVsNeutrinoKinematics] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();

    ROOT::RDF::RNode node_mc = filter_by_sample_mask(rdf, mask_mc);
    node_mc = filter_not_sample_mask(node_mc, mask_ext);

    const std::vector<VarSpec> vars = {
        {"nu_E",            30, 0.0, 10.0, "True #nu energy E_{#nu} [GeV]"},
        {"Q2",              30, 0.0, 5.0,  "True Q^{2} [GeV^{2}]"},
        {"kin_W",           30, 0.0, 5.0,  "True W [GeV]"},
        {"bjorken_x",       30, 0.0, 1.0,  "True Bjorken-x"},
        {"inelasticity_y",  30, 0.0, 1.0,  "True inelasticity y"},
        {"nu_theta",        30, 0.0, 3.2,  "True #nu #theta [rad]"},
        {"nu_pt",           30, 0.0, 5.0,  "True #nu p_{T} [GeV/c]"},
    };

    int rc = 0;
    for (const auto &v : vars)
    {
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
        opt.out_dir.clear();
        opt.image_format.clear();

        EfficiencyPlot eff(spec, opt, cfg);
        const int one = eff.compute(node_mc, denom_sel, pass_sel, extra_sel);
        if (one == 0)
        {
            const std::string tag = Plotter::sanitise(v.expr);
            rc = std::max(rc, eff.draw_and_save("pr_eff_vs_" + tag));
        }
        else
        {
            rc = std::max(rc, one);
        }
    }

    return rc;
}
