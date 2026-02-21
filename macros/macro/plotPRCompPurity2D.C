// plot/macro/plotPRCompPurity2D.C
//
// Plot 2D histograms of pattern-recognition completeness vs purity for:
//   - pion   (pr_pi_completeness vs pr_pi_purity)
//   - proton (pr_p_completeness  vs pr_p_purity)
//   - muon   (pr_mu_completeness vs pr_mu_purity)
//
// Default selections mirror the PR-efficiency macros:
//   - Base denom (eligible truth events):
//       is_nu_mu_cc && nu_vtx_in_fv && (lam_pdg==3122) && (p_p>0) && (pi_p>0)
//   - Default requires a valid assignment for meaningful comp/pur:
//       pr_valid_assignment
//
// Run with:
//   ./heron macro plotPRCompPurity2D.C
//   ./heron macro plotPRCompPurity2D.C 'plotPRCompPurity2D("./scratch/out/event_list_myana.root","sel_triggered_slice")'
//   ./heron macro plotPRCompPurity2D.C 'plotPRCompPurity2D("./scratch/out/event_list_myana.root","sel_triggered_slice","pr_valid_assignment")'
//   ./heron macro plotPRCompPurity2D.C 'plotPRCompPurity2D("./scratch/out/event_list_myana.root","sel_triggered_slice","true")'  // include even unassigned (still filters finite)
//
// Output:
//   Saves:
//     - one combined 1x3 canvas:  pr_comp_pur_2d_all.<fmt>
//     - and one per particle:     pr_comp_pur_2d_{mu,p,pi}.<fmt>
//   to:
//     $HERON_PLOT_DIR (default: ./scratch/plots)
//   with format:
//     $HERON_PLOT_FORMAT (default: pdf)

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventListIO.hh"
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

struct PartSpec
{
    std::string tag;      // used in filenames
    std::string label;    // plot label
    std::string comp;     // x expr
    std::string pur;      // y expr
};

struct PlotCfg
{
    int nbins = 60;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;

    // Optional “paper cut” guides (same as default PR-eff macro thresholds)
    double cut_comp = 0.1;
    double cut_pur = 0.5;

    // Require values to be in [0,1]? (filters out sentinel values like -999)
    bool filter_unit_interval = true;

    // Make low-occupancy structure visible.
    bool logz = true;
    double logz_min = 0.5; // must be > 0 for log scale
};

std::string finite_pair_sel(const std::string &x, const std::string &y)
{
    // NaN check: x==x is false iff x is NaN
    return "(" + x + " == " + x + ") && (" + y + " == " + y + ")";
}

std::string unit_interval_sel(const std::string &x, const std::string &y,
                              double xmin, double xmax, double ymin, double ymax)
{
    std::ostringstream os;
    os << "(" << x << " >= " << xmin << " && " << x << " <= " << xmax << ")"
       << " && "
       << "(" << y << " >= " << ymin << " && " << y << " <= " << ymax << ")";
    return os.str();
}

void print_comp_pur_stats(TH2D &h,
                          const std::string &label,
                          const PlotCfg &cfg)
{
    const int nx = h.GetNbinsX();
    const int ny = h.GetNbinsY();

    // Strict ">" behaviour even if cut sits on a bin edge.
    const int bx = h.GetXaxis()->FindBin(std::nextafter(cfg.cut_comp, cfg.xmax));
    const int by = h.GetYaxis()->FindBin(std::nextafter(cfg.cut_pur, cfg.ymax));

    const double n_tot = h.Integral(1, nx, 1, ny);
    const double n_comp = h.Integral(bx, nx, 1, ny);
    const double n_pur = h.Integral(1, nx, by, ny);
    const double n_both = h.Integral(bx, nx, by, ny);

    auto pct = [](double num, double den) -> double {
        return (den > 0.0) ? (100.0 * num / den) : 0.0;
    };

    std::ostringstream os;
    os << "[plotPRCompPurity2D] " << label << " comp/pur stats\n";
    os << "  entries(fills): " << static_cast<long long>(h.GetEntries())
       << ", sumw(in-range): " << n_tot << "\n";

    os << std::fixed << std::setprecision(1);
    os << "  comp > " << cfg.cut_comp << ": " << pct(n_comp, n_tot) << "% (sumw=" << n_comp << ")\n";
    os << "  pur  > " << cfg.cut_pur << ": " << pct(n_pur, n_tot) << "% (sumw=" << n_pur << ")\n";
    os << "  both: " << pct(n_both, n_tot) << "% (sumw=" << n_both << ")\n";
    std::cout << os.str();
}

void draw_comp_pur_hist(TH2D &h,
                        const PlotCfg &cfg)
{
    // Axis labels only (no title text block)
    h.SetTitle(";Completeness;Purity");
    h.GetXaxis()->SetRangeUser(cfg.xmin, cfg.xmax);
    h.GetYaxis()->SetRangeUser(cfg.ymin, cfg.ymax);
    h.GetZaxis()->SetTitle("Events");

    // Make the color bar fit
    if (gPad)
    {
        gPad->SetRightMargin(0.14);
        gPad->SetLogz(cfg.logz ? 1 : 0);
    }

    if (cfg.logz)
        h.SetMinimum(cfg.logz_min);

    h.Draw("COLZ");

    // Draw cut guide lines (heap-allocated so ROOT can own them via pad primitives)
    {
        TLine *lx = new TLine(cfg.cut_comp, cfg.ymin, cfg.cut_comp, cfg.ymax);
        lx->SetLineStyle(2);
        lx->SetLineWidth(3);
        lx->SetLineColor(kRed);
        lx->Draw("same");

        TLine *ly = new TLine(cfg.xmin, cfg.cut_pur, cfg.xmax, cfg.cut_pur);
        ly->SetLineStyle(2);
        ly->SetLineWidth(3);
        ly->SetLineColor(kRed);
        ly->Draw("same");
    }
}

} // namespace

int plotPRCompPurity2D(const std::string &samples_tsv = "",
                       const std::string &extra_sel = "true",
                       const std::string &assignment_sel = "pr_valid_assignment",
                       const std::string &denom_sel =
                           "is_nu_mu_cc"
                           " && nu_vtx_in_fv"
                           " && (lam_pdg==3122)"
                           " && (p_p>0.0)"
                           " && (pi_p>0.0)",
                       int nbins = 60,
                       double cut_comp = 0.1,
                       double cut_pur = 0.5,
                       bool filter_unit_interval = true,
                       bool logz = true)
{
  return heron::macro::run_with_guard("plotPRCompPurity2D", [&]() -> int {
    ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotPRCompPurity2D] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotPRCompPurity2D] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();

    auto filter_by_mask = [](ROOT::RDF::RNode n,
                             std::shared_ptr<const std::vector<char>> mask) -> ROOT::RDF::RNode {
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

    // Apply the same kind of “base” selections as the PR-eff macros.
    ROOT::RDF::RNode base = node_mc;
    if (!extra_sel.empty())
        base = base.Filter(extra_sel);
    if (!denom_sel.empty())
        base = base.Filter(denom_sel);
    if (!assignment_sel.empty())
        base = base.Filter(assignment_sel);

    const std::string out_dir = getenv_or("HERON_PLOT_DIR", "./scratch/plots");
    const std::string out_fmt = getenv_or("HERON_PLOT_FORMAT", "pdf");
    gSystem->mkdir(out_dir.c_str(), /*recursive=*/true);

    // Configure plot appearance
    PlotCfg cfg;
    cfg.nbins = nbins;
    cfg.cut_comp = cut_comp;
    cfg.cut_pur = cut_pur;
    cfg.filter_unit_interval = filter_unit_interval;
    cfg.logz = logz;

    const std::vector<PartSpec> parts = {
        {"mu", "Muon",   "pr_mu_completeness", "pr_mu_purity"},
        {"p",  "Proton", "pr_p_completeness",  "pr_p_purity"},
        {"pi", "Pion",   "pr_pi_completeness", "pr_pi_purity"},
    };

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);

    // Build histograms (keep RResultPtr alive for combined canvas)
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h2;
    h2.reserve(parts.size());

    for (const auto &p : parts)
    {
        ROOT::RDF::RNode n = base;

        // Finite checks
        n = n.Filter(finite_pair_sel(p.comp, p.pur));

        // Optional unit-interval enforcement (filters sentinel values)
        if (cfg.filter_unit_interval)
            n = n.Filter(unit_interval_sel(p.comp, p.pur, cfg.xmin, cfg.xmax, cfg.ymin, cfg.ymax));

        const std::string hname = "h2_" + sanitize_for_filename(p.tag);
        auto hist = n.Histo2D({hname.c_str(), "", cfg.nbins, cfg.xmin, cfg.xmax, cfg.nbins, cfg.ymin, cfg.ymax},
                             p.comp, p.pur);
        h2.push_back(hist);

        // Print stats once per particle to terminal (no on-plot box)
        print_comp_pur_stats(*hist, p.label, cfg);

        // Save individual plot
        TCanvas c(("c2_" + p.tag).c_str(), "", 900, 800);
        draw_comp_pur_hist(*hist, cfg);

        const std::string out_path = out_dir + "/pr_comp_pur_2d_" + p.tag + "." + out_fmt;
        c.SaveAs(out_path.c_str());
        std::cout << "[plotPRCompPurity2D] wrote " << out_path << "\n";
    }

    // Save combined 1x3 canvas
    {
        TCanvas c_all("c2_all", "", 1500, 520);
        c_all.Divide(3, 1);

        for (size_t i = 0; i < parts.size(); ++i)
        {
            c_all.cd(static_cast<int>(i) + 1);
            draw_comp_pur_hist(*h2[i], cfg);
        }

        const std::string out_path = out_dir + "/pr_comp_pur_2d_all." + out_fmt;
        c_all.SaveAs(out_path.c_str());
        std::cout << "[plotPRCompPurity2D] wrote " << out_path << "\n";
    }

    return 0;

  });
}
