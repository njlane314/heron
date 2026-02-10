// plot/macro/plotSignalCoverageTruthKinematics.C
//
// Make the "truth coverage" distributions shown in the draft screenshot:
//
//   Panel A: Truth E_{#nu} distribution
//   Panel B: Truth W distribution
//   Panel D (optional): 2D truth-kinematics coverage plots
//
// The 1D plots are stacked by analysis channel via Plotter (same machinery used in
// plotStackedHistTrueVertex.C), with signal overlay enabled.
//
// Run with:
//   ./nuxsec macro plotSignalCoverageTruthKinematics.C
//   ./nuxsec macro plotSignalCoverageTruthKinematics.C \
//        'plotSignalCoverageTruthKinematics("./scratch/out/event_list_myana.root","sel_cc0pi")'
//
// Notes:
//   - This macro expects an event-list ROOT file (event_list_<analysis>.root).
//   - Truth variables are MC-only; by default this macro plots MC only.
//   - Binning is fixed ("common binning") via the VarSpec tables below.
//
// Output:
//   - 1D stacked plots follow Plotter/PlotEnv defaults (NUXSEC_PLOT_DIR / NUXSEC_PLOT_FORMAT).
//   - 2D plots are written explicitly to $NUXSEC_PLOT_DIR (default: ./scratch/plots).
//
// You may pass an explicit selection expression via `extra_sel`.

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
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

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

bool implicit_mt_enabled()
{
    const char *env = std::getenv("NUXSEC_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}

struct Var1D
{
    // `name` controls output naming inside Plotter (hist name / file stem).
    std::string name;
    // `expr` is the RDataFrame column/expression to plot.
    std::string expr;
    int nbins = 30;
    double xmin = 0.0;
    double xmax = 1.0;
    std::string x_title;
};

struct Var2D
{
    std::string name;
    std::string x_expr;
    std::string y_expr;
    int nx = 30;
    double xmin = 0.0;
    double xmax = 1.0;
    int ny = 30;
    double ymin = 0.0;
    double ymax = 1.0;
    std::string x_title;
    std::string y_title;
};

void draw_truth_2d(ROOT::RDF::RNode node,
                   const Var2D &v,
                   const std::string &weight,
                   const std::string &out_dir,
                   const std::string &out_fmt)
{
    // Drop NaNs explicitly for both axes.
    ROOT::RDF::RNode n2 = node.Filter(v.x_expr + " == " + v.x_expr)
                             .Filter(v.y_expr + " == " + v.y_expr);

    const auto n_evt = n2.Count().GetValue();
    if (n_evt == 0)
    {
        std::cout << "[plotSignalCoverageTruthKinematics] skip 2D " << v.name
                  << " (no entries after selections)\n";
        return;
    }

    const std::string hname = "h2_" + sanitize_for_filename(v.name);
    auto h2 = n2.Histo2D({hname.c_str(), "", v.nx, v.xmin, v.xmax, v.ny, v.ymin, v.ymax},
                         v.x_expr, v.y_expr, weight);

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);

    TCanvas c(("c2_" + sanitize_for_filename(v.name)).c_str(), "", 900, 750);
    c.SetRightMargin(0.14);
    c.SetGrid();

    const std::string title = ";" + (v.x_title.empty() ? v.x_expr : v.x_title) +
                              ";" + (v.y_title.empty() ? v.y_expr : v.y_title) +
                              ";Events";
    h2->SetTitle(title.c_str());
    h2->Draw("COLZ");

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.032);
    std::ostringstream os;
    os << "Entries: " << n_evt;
    lat.DrawLatex(0.14, 0.92, os.str().c_str());

    gSystem->mkdir(out_dir.c_str(), /*recursive=*/true);

    const std::string out_path = out_dir + "/" + sanitize_for_filename(v.name) + "." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plotSignalCoverageTruthKinematics] wrote " << out_path << "\n";
}

} // namespace

int plotSignalCoverageTruthKinematics(const std::string &samples_tsv = "",
                                      // Optional additional selection expression (e.g. "sel_cc0pi").
                                      const std::string &extra_sel = "true",
                                      bool make_2d = true)
{
    if (implicit_mt_enabled())
    {
        ROOT::EnableImplicitMT();
        std::cout << "[plotSignalCoverageTruthKinematics] ROOT implicit MT enabled (NUXSEC_PLOT_IMT != 0)\n";
    }

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotSignalCoverageTruthKinematics] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotSignalCoverageTruthKinematics] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    // --- Select MC only (exclude EXT) ---
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

    // --- Build a single MC Entry for Plotter ---
    std::vector<Entry> entries;
    entries.reserve(1);

    std::vector<const Entry *> mc;

    ProcessorEntry rec_mc;
    rec_mc.source = Type::kMC;

    entries.emplace_back(make_entry(std::move(node_mc), rec_mc));
    Entry &e_mc = entries.back();
    mc.push_back(&e_mc);

    if (!extra_sel.empty())
    {
        std::cout << "[plotSignalCoverageTruthKinematics] selection=" << extra_sel << "\n";
        e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(extra_sel);
    }

    // --- Configure Plotter options (match existing stacked-hist style) ---
    Plotter plotter;
    auto &opt = plotter.options();
    opt.use_log_y = false;
    opt.legend_on_top = true;
    opt.annotate_numbers = true;
    opt.overlay_signal = true;
    opt.show_ratio = false;
    opt.show_ratio_band = false;
    opt.adaptive_binning = false;          // fixed binning (paper wants common binning)
    opt.adaptive_fold_overflow = true;
    opt.signal_channels = Channels::signal_keys();
    opt.y_title = "Events/bin";
    opt.run_numbers = {"1"};
    opt.image_format = getenv_or("NUXSEC_PLOT_FORMAT", "pdf");

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();

    const std::string weight = "w_nominal";

    // --- Panel A/B: 1D truth distributions (stacked by analysis channel) ---
    const std::vector<Var1D> vars_1d = {
        // Panel A
        {"truth_nu_E", "nu_E",   30, 0.0, 10.0, "True #nu energy E_{#nu} [GeV]"},
        // Panel B
        {"truth_W",    "kin_W",  30, 0.0, 5.0,  "True W [GeV]"},
    };

    for (const auto &v : vars_1d)
    {
        opt.x_title = v.x_title.empty() ? v.expr : v.x_title;

        // `make_spec(name, ...)` controls plot/hist naming; `spec.expr` controls what is drawn.
        TH1DModel spec = make_spec(v.name, v.nbins, v.xmin, v.xmax, weight);
        spec.expr = v.expr;
        spec.sel = Preset::Empty;

        plotter.draw_stack(spec, mc);
    }

    // --- Optional: Panel D-style 2D truth coverage plots (signal-only recommended) ---
    if (make_2d)
    {
        // "Truth signal" definition used for the efficiency plots (adjust if needed).
        const std::string truth_signal_sel =
            "is_nu_mu_cc"
            " && nu_vtx_in_fv"
            " && (lam_pdg==3122)"
            " && (p_p>0.0)"
            " && (pi_p>0.0)";

        ROOT::RDF::RNode node_sig = e_mc.selection.nominal.node.Filter(truth_signal_sel);

        const std::string out_dir = getenv_or("NUXSEC_PLOT_DIR", "./scratch/plots");
        const std::string out_fmt = getenv_or("NUXSEC_PLOT_FORMAT", "pdf");

        const std::vector<Var2D> vars_2d = {
            {"truth2d_W_vs_nu_E",     "nu_E",      "kin_W",
                                     40, 0.0, 10.0, 40, 0.0, 5.0,
                                     "True E_{#nu} [GeV]", "True W [GeV]"},
            {"truth2d_x_vs_W",        "kin_W",     "bjorken_x",
                                     40, 0.0, 5.0,  40, 0.0, 1.0,
                                     "True W [GeV]", "True Bjorken-x"},
            {"truth2d_Q2_vs_nu_E",    "nu_E",      "Q2",
                                     40, 0.0, 10.0, 40, 0.0, 5.0,
                                     "True E_{#nu} [GeV]", "True Q^{2} [GeV^{2}]"},
        };

        for (const auto &v : vars_2d)
        {
            draw_truth_2d(node_sig, v, weight, out_dir, out_fmt);
        }
    }

    return 0;
}
