// plot/macro/plotPREffVsLambdaKinematics.C
//
// Plot pattern-recognition "assignment efficiency" vs Lambda kinematics.
//
// Efficiency definition (default):
//   denom_sel: mu_true_hits>0 && p_true_hits>0 && pi_true_hits>0
//   pass_sel:  pr_valid_assignment
//
// Run with:
//   ./nuxsec macro plotPREffVsLambdaKinematics.C
//   ./nuxsec macro plotPREffVsLambdaKinematics.C 'plotPREffVsLambdaKinematics("./scratch/out/event_list_myana.root","sel_triggered_slice")'
//   ./nuxsec macro plotPREffVsLambdaKinematics.C 'plotPREffVsLambdaKinematics("./scratch/out/event_list_myana.root","sel_triggered_slice","pr_valid_assignment && pr_mu_completeness>0.5 && pr_p_completeness>0.5 && pr_pi_completeness>0.5")'
//
// Output:
//   Saves one plot per x-variable to:
//     $NUXSEC_PLOT_DIR (default: ./scratch/plots)
//   with format:
//     $NUXSEC_PLOT_FORMAT (default: pdf)

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventListIO.hh"
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

std::string getenv_or(const char *key, const std::string &fallback)
{
    const char *v = std::getenv(key);
    if (!v || std::string(v).empty())
        return fallback;
    return std::string(v);
}

std::string sanitize_for_filename(const std::string &s)
{
    std::string out;
    out.reserve(s.size());
    for (const char c : s)
    {
        const bool ok =
            (c >= 'a' && c <= 'z') ||
            (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') ||
            (c == '_');
        out.push_back(ok ? c : '_');
    }
    return out;
}

struct Axis
{
    int nbins = 30;
    double xmin = 0.0;
    double xmax = 1.0;
};

Axis build_axis(ROOT::RDF::RNode node_finite,
                const std::string &expr,
                int fallback_nbins,
                double fallback_xmin,
                double fallback_xmax)
{
    Axis out;
    out.nbins = fallback_nbins;
    out.xmin = fallback_xmin;
    out.xmax = fallback_xmax;

    const auto n = node_finite.Count().GetValue();
    if (n == 0)
        return out;

    const double vmin = static_cast<double>(node_finite.Min(expr).GetValue());
    const double vmax = static_cast<double>(node_finite.Max(expr).GetValue());

    if (!std::isfinite(vmin) || !std::isfinite(vmax) || vmax <= vmin)
        return out;

    const double span = vmax - vmin;
    const double pad = (span > 0.0 ? 0.05 * span : 1.0);
    out.xmin = vmin - pad;
    out.xmax = vmax + pad;
    return out;
}

struct VarSpec
{
    std::string expr;
    int nbins = 30;
    double xmin = 0.0;
    double xmax = 1.0;
    std::string x_title;
};

int draw_efficiency(ROOT::RDF::RNode mc_base,
                    const VarSpec &v,
                    const std::string &extra_sel,
                    const std::string &denom_sel,
                    const std::string &pass_sel,
                    const std::string &out_dir,
                    const std::string &out_fmt)
{
    // Build denom and numerator nodes.
    ROOT::RDF::RNode denom = mc_base;
    if (!extra_sel.empty())
        denom = denom.Filter(extra_sel);
    if (!denom_sel.empty())
        denom = denom.Filter(denom_sel);

    // Exclude NaNs for this variable explicitly.
    ROOT::RDF::RNode denom_finite = denom.Filter(v.expr + " == " + v.expr);

    const auto n_denom = denom_finite.Count().GetValue();
    if (n_denom == 0)
    {
        std::cout << "[plotPREffVsLambdaKinematics] skip " << v.expr
                  << " (no denom entries after selection)\n";
        return 0;
    }

    ROOT::RDF::RNode numer = denom_finite;
    if (!pass_sel.empty())
        numer = numer.Filter(pass_sel);

    const auto n_numer = numer.Count().GetValue();

    const Axis ax = build_axis(denom_finite, v.expr, v.nbins, v.xmin, v.xmax);

    const std::string tag = sanitize_for_filename(v.expr);
    const std::string hden_name = "h_den_" + tag;
    const std::string hnum_name = "h_num_" + tag;

    auto h_den = denom_finite.Histo1D({hden_name.c_str(), "", ax.nbins, ax.xmin, ax.xmax}, v.expr);
    auto h_num = numer.Histo1D({hnum_name.c_str(), "", ax.nbins, ax.xmin, ax.xmax}, v.expr);

    if (!TEfficiency::CheckConsistency(*h_num, *h_den))
    {
        std::cerr << "[plotPREffVsLambdaKinematics] TEfficiency consistency check failed for "
                  << v.expr << "\n";
        return 1;
    }

    TEfficiency eff(*h_num, *h_den);
    eff.SetStatisticOption(TEfficiency::kFCP); // Clopper-Pearson

    std::unique_ptr<TGraphAsymmErrors> g(eff.CreateGraph());
    if (!g)
    {
        std::cerr << "[plotPREffVsLambdaKinematics] failed to create graph for " << v.expr << "\n";
        return 1;
    }

    const std::string title = ";" + (v.x_title.empty() ? v.expr : v.x_title) + ";Pattern recognition efficiency";
    g->SetTitle(title.c_str());
    g->GetYaxis()->SetRangeUser(0.0, 1.05);

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);

    TCanvas c(("c_" + tag).c_str(), "", 900, 700);
    c.SetGrid();

    g->Draw("AP");

    // Annotate counts + overall efficiency.
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.032);

    std::ostringstream os1;
    os1 << "Denom: " << n_denom << "   Numer: " << n_numer;
    lat.DrawLatex(0.14, 0.92, os1.str().c_str());

    std::ostringstream os2;
    const double eff_all = (n_denom > 0 ? static_cast<double>(n_numer) / static_cast<double>(n_denom) : 0.0);
    os2 << "Overall eff: " << eff_all;
    lat.DrawLatex(0.14, 0.88, os2.str().c_str());

    gSystem->mkdir(out_dir.c_str(), /*recursive=*/true);

    const std::string out_path = out_dir + "/pr_eff_vs_" + tag + "." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plotPREffVsLambdaKinematics] wrote " << out_path << "\n";

    return 0;
}

} // namespace

int plotPREffVsLambdaKinematics(const std::string &samples_tsv = "",
                                const std::string &extra_sel = "true",
                                const std::string &pass_sel = "pr_valid_assignment",
                                const std::string &denom_sel = "mu_true_hits>0 && p_true_hits>0 && pi_true_hits>0")
{
    ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotPREffVsLambdaKinematics] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotPREffVsLambdaKinematics] input is not an event list root file: " << list_path << "\n";
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

    const std::string out_dir = getenv_or("NUXSEC_PLOT_DIR", "./scratch/plots");
    const std::string out_fmt = getenv_or("NUXSEC_PLOT_FORMAT", "pdf");

    // Lambda + decay + muon kinematics available in LambdaAnalysis_tool.
    const std::vector<VarSpec> vars = {
        {"lam_p_mag",          30, 0.0, 3.0,  "#Lambda |p| [GeV/c]"},
        {"lam_E",              30, 0.0, 5.0,  "#Lambda E [GeV]"},
        {"lam_ct",             30, 0.0, 100.0,"#Lambda c#tau proxy [cm]"},
        {"lam_decay_sep",      30, 0.0, 200.0,"#Lambda decay sep. to #nu vtx [cm]"},
        {"p_p",                30, 0.0, 3.0,  "p |p| [GeV/c]"},
        {"pi_p",               30, 0.0, 3.0,  "#pi |p| [GeV/c]"},
        {"ppi_opening_angle",  30, 0.0, 3.2,  "Opening angle(p,#pi) [rad]"},
        {"mu_p",               30, 0.0, 5.0,  "#mu |p| [GeV/c]"},
        {"mu_theta",           30, 0.0, 3.2,  "#theta(#nu,#mu) [rad]"},
    };

    int rc = 0;
    for (const auto &v : vars)
    {
        const int one = draw_efficiency(node_mc, v, extra_sel, denom_sel, pass_sel, out_dir, out_fmt);
        rc = std::max(rc, one);
    }

    return rc;
}
