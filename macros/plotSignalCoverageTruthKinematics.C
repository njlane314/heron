// plot/macro/plotSignalCoverageTruthKinematics.C
//
// Make the "truth coverage" distributions shown in the draft screenshot:
//
//   Panel A: Truth E_{#nu} distribution
//   Panel B: Truth W distribution
//   Panel D (optional): 2D truth-kinematics coverage plots
//
// The 1D plots are stacked by analysis channel via Plotter, with signal overlay
// enabled.
//
// Run with:
//   ./heron macro plotSignalCoverageTruthKinematics.C
//   ./heron macro plotSignalCoverageTruthKinematics.C \
//        'plotSignalCoverageTruthKinematics("./scratch/out/event_list_myana.root","sel_cc0pi")'
//   ./heron macro plotSignalCoverageTruthKinematics.C \
//        'plotSignalCoverageTruthKinematics("./scratch/out/event_list_myana.root",
//                                          "sel_cc0pi",
//                                          true,  /*make_2d*/
//                                          true)' /*signal_only_2d*/
//
// Notes:
//   - This macro expects an event-list ROOT file (event_list_<analysis>.root).
//   - Truth variables are MC-only; by default this macro plots MC only.
//   - Binning is fixed ("common binning") via the VarSpec tables below.
//
// Output:
//   - 1D stacked plots follow Plotter/PlotEnv defaults (HERON_PLOT_DIR / HERON_PLOT_FORMAT).
//   - 2D plots are written explicitly to $HERON_PLOT_DIR (default: ./scratch/plots).
//
// You may pass an explicit selection expression via `extra_sel`.

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
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "include/MacroGuard.hh"
#include "include/MacroEnv.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{



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
    const char *env = std::getenv("HERON_PLOT_IMT");
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
    std::string y_title;
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

std::string finite_pair_sel(const std::string &x, const std::string &y)
{
    // NaN check: x==x is false iff x is NaN
    return "(" + x + " == " + x + ") && (" + y + " == " + y + ")";
}

std::string in_range_sel(const Var2D &v)
{
    // Explicitly keep only values that will land in visible bins (avoid silent under/overflow-only fills).
    std::ostringstream os;
    os << "(" << v.x_expr << " >= " << v.xmin << " && " << v.x_expr << " <= " << v.xmax << ")"
       << " && "
       << "(" << v.y_expr << " >= " << v.ymin << " && " << v.y_expr << " <= " << v.ymax << ")";
    return os.str();
}

double min_positive_bin_content(const TH2D &h)
{
    const int nx = h.GetNbinsX();
    const int ny = h.GetNbinsY();
    double min_pos = std::numeric_limits<double>::infinity();
    for (int ix = 1; ix <= nx; ++ix)
        for (int iy = 1; iy <= ny; ++iy)
        {
            const double z = h.GetBinContent(ix, iy);
            if (std::isfinite(z) && z > 0.0 && z < min_pos)
                min_pos = z;
        }
    return min_pos;
}

bool any_negative_bin(const TH2D &h)
{
    const int nx = h.GetNbinsX();
    const int ny = h.GetNbinsY();
    for (int ix = 1; ix <= nx; ++ix)
        for (int iy = 1; iy <= ny; ++iy)
            if (h.GetBinContent(ix, iy) < 0.0)
                return true;
    return false;
}

void draw_truth_2d(ROOT::RDF::RNode node,
                   const Var2D &v,
                   const std::string &weight,
                   const std::string &out_dir,
                   const std::string &out_fmt)
{
    const bool use_logz = true;
    const double logz_min_floor = 1e-3;

    // Drop NaNs and enforce plotted range (avoid filling only under/overflow).
    ROOT::RDF::RNode n2 = node.Filter(finite_pair_sel(v.x_expr, v.y_expr))
                             .Filter(in_range_sel(v));

    // If weighting, require finite weights; if logz, require positive weights.
    // (Negative weights + logz will produce misleading/blank plots.)
    if (!weight.empty() && weight != "1" && weight != "1.0")
    {
        n2 = n2.Filter(weight + " == " + weight);
        if (use_logz)
            n2 = n2.Filter(weight + " > 0");
    }

    const std::string hname = "h2_" + sanitize_for_filename(v.name);
    ROOT::RDF::RResultPtr<TH2D> h2;
    if (weight.empty() || weight == "1" || weight == "1.0")
    {
        h2 = n2.Histo2D({hname.c_str(), "", v.nx, v.xmin, v.xmax, v.ny, v.ymin, v.ymax},
                        v.x_expr, v.y_expr);
    }
    else
    {
        h2 = n2.Histo2D({hname.c_str(), "", v.nx, v.xmin, v.xmax, v.ny, v.ymin, v.ymax},
                        v.x_expr, v.y_expr, weight);
    }

    const auto n_entries = static_cast<long long>(h2->GetEntries());
    if (n_entries == 0)
    {
        std::cout << "[plotSignalCoverageTruthKinematics] skip 2D " << v.name
                  << " (no entries after selections)\n";
        return;
    }

    // Print stats to terminal (match plotPRCompPurity2D style; no on-plot box)
    {
        const int nx = h2->GetNbinsX();
        const int ny = h2->GetNbinsY();
        const double sumw = h2->Integral(1, nx, 1, ny);
        std::cout << "[plotSignalCoverageTruthKinematics] " << v.name << " 2D stats\n"
                  << "  entries(fills): " << n_entries << ", sumw(in-range): " << sumw << "\n";
        if (sumw <= 0.0)
        {
            std::cout << "[plotSignalCoverageTruthKinematics] skip 2D " << v.name
                      << " (no in-range content)\n";
            return;
        }
    }

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(255);
    gStyle->SetPalette();

    TCanvas c(("c2_" + sanitize_for_filename(v.name)).c_str(), "", 920, 800);
    c.cd(); // ensure gPad is set

    // Axis labels only (no title text block); z-title set explicitly
    const std::string xlab = (v.x_title.empty() ? v.x_expr : v.x_title);
    const std::string ylab = (v.y_title.empty() ? v.y_expr : v.y_title);
    h2->SetTitle((";" + xlab + ";" + ylab).c_str());
    h2->GetZaxis()->SetTitle("Events");
    h2->GetXaxis()->SetRangeUser(v.xmin, v.xmax);
    h2->GetYaxis()->SetRangeUser(v.ymin, v.ymax);
    h2->GetXaxis()->SetTitleOffset(1.15);
    h2->GetYaxis()->SetTitleOffset(1.35);
    h2->GetZaxis()->SetTitleOffset(1.2);
    h2->GetXaxis()->SetNdivisions(505);
    h2->GetYaxis()->SetNdivisions(505);

    // Make the color bar fit; apply logz at the pad level (like plotPRCompPurity2D)
    const bool has_neg = any_negative_bin(*h2);
    const bool do_logz = use_logz && !has_neg;
    if (gPad)
    {
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.12);
        gPad->SetRightMargin(0.18);
        gPad->SetTopMargin(0.08);
        gPad->SetLogz(do_logz ? 1 : 0);
    }
    if (do_logz)
    {
        const double min_pos = min_positive_bin_content(*h2);
        const double zmin = (std::isfinite(min_pos) ? std::max(logz_min_floor, 0.5 * min_pos) : logz_min_floor);
        h2->SetMinimum(zmin);
    }

    h2->Draw("COLZ");

    gSystem->mkdir(out_dir.c_str(), /*recursive=*/true);

    const std::string out_path = out_dir + "/" + sanitize_for_filename(v.name) + "." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plotSignalCoverageTruthKinematics] wrote " << out_path << "\n";
}

} // namespace

int plotSignalCoverageTruthKinematics(const std::string &samples_tsv = "",
                                      // Optional additional selection expression (e.g. "sel_cc0pi").
                                      const std::string &extra_sel = "true",
                                      bool make_2d = true,
                                      bool signal_only_2d = false)
{
  return heron::macro::run_with_guard("plotSignalCoverageTruthKinematics", [&]() -> int {
    if (implicit_mt_enabled())
    {
        ROOT::EnableImplicitMT();
        std::cout << "[plotSignalCoverageTruthKinematics] ROOT implicit MT enabled (HERON_PLOT_IMT != 0)\n";
    }

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotSignalCoverageTruthKinematics] input=" << list_path << "\n";

    if (!heron::macro::looks_like_event_list_root(list_path))
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
    opt.signal_channels = Channels::signal_keys();
    opt.y_title = "Events / bin";
    opt.run_numbers = {"1"};
    opt.image_format = heron::macro::getenv_or("HERON_PLOT_FORMAT", "pdf");

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();

    const std::string weight = "w_nominal";

    // --- Panel A/B: 1D truth distributions (stacked by analysis channel) ---
    const std::vector<Var1D> vars_1d = {
        // Panel A
        {"truth_nu_E", "nu_E",   30, 0.0, 10.0,
         "True #nu energy E_{#nu} [GeV]", "Events / #DeltaE_{#nu} [GeV]"},
        // Panel B
        {"truth_W",    "kin_W",  30, 0.0, 5.0,
         "True W [GeV]", "Events / #DeltaW [GeV]"},
    };

    for (const auto &v : vars_1d)
    {
        opt.x_title = v.x_title.empty() ? v.expr : v.x_title;
        opt.y_title = v.y_title.empty() ? "Events / bin" : v.y_title;

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

        // By default, plot *all selected events* (like plotPRCompPurity2D does after its base selections).
        // If you want signal-only 2D, pass signal_only_2d=true.
        ROOT::RDF::RNode node_2d = e_mc.selection.nominal.node;
        if (signal_only_2d)
        {
            node_2d = node_2d.Filter(truth_signal_sel);
            std::cout << "[plotSignalCoverageTruthKinematics] 2D: truth-signal-only selection enabled\n";
        }
        else
        {
            std::cout << "[plotSignalCoverageTruthKinematics] 2D: plotting all selected events (no truth-signal filter)\n";
        }

        const std::string out_dir = heron::macro::getenv_or("HERON_PLOT_DIR", "./scratch/plots");
        const std::string out_fmt = heron::macro::getenv_or("HERON_PLOT_FORMAT", "pdf");

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
            draw_truth_2d(node_2d, v, weight, out_dir, out_fmt);
        }
    }

    return 0;

  });
}
