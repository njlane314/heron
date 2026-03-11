// macros/plot_nominal_weight_signal_bkg_unstacked.C
//
// Overlay (unstacked) histograms of the nominal event weight for:
//   - signal (MC only)
//   - total background (MC background + EXT by default)
// after a base selection (default: sel_muon).
//
// Important:
//   The x-axis is w_nominal itself, so the histogram is filled UNWEIGHTED.
//   In other words, w_nominal is treated as the observable to plot, not as
//   the histogram weight.
//
// Run with:
//   ./heron macro plot_nominal_weight_signal_bkg_unstacked.C
//   ./heron macro plot_nominal_weight_signal_bkg_unstacked.C \
//     'plot_nominal_weight_signal_bkg_unstacked("./scratch/out/event_list_myana.root")'
//
// Example with explicit range:
//   ./heron macro plot_nominal_weight_signal_bkg_unstacked.C \
//     'plot_nominal_weight_signal_bkg_unstacked("", "sel_muon", "is_signal", "w_nominal", 80, 0.0, 5.0, true, true, "nominal_weight_after_sel_muon")'

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

void set_hist_style(TH1D &h, int color, int marker)
{
    h.SetLineColor(color);
    h.SetMarkerColor(color);
    h.SetLineWidth(3);
    h.SetMarkerStyle(marker);
    h.SetMarkerSize(0.9);
    h.SetFillStyle(0);
}

} // namespace

int plot_nominal_weight_signal_bkg_unstacked(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &signal_sel = "is_signal",
    const std::string &weight_col = "w_nominal",
    int nbins = 80,
    double xmin = -1.0,
    double xmax = -1.0,
    bool include_ext_in_background = true,
    bool use_logy = true,
    const std::string &output_stem = "nominal_weight_signal_bkg_unstacked")
{
    return heron::macro::run_with_guard("plot_nominal_weight_signal_bkg_unstacked", [&]() -> int {
        ROOT::EnableImplicitMT();
        TH1::SetDefaultSumw2();

        const std::string input_path =
            event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_nominal_weight_signal_bkg_unstacked] input="
                  << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_nominal_weight_signal_bkg_unstacked] input is not an event-list ROOT file: "
                      << input_path << "\n";
            return 1;
        }

        if (nbins < 1)
            nbins = 1;

        EventListIO el(input_path);

        ROOT::RDF::RNode base = SelectionService::decorate(el.rdf())
            // Treat w_nominal as the observable to plot.
            // Do NOT use it as the histogram weight.
            .Define("__plot_x__", weight_col)
            .Filter([](double x) { return std::isfinite(x); }, {"__plot_x__"});

        if (!base_sel.empty())
        {
            if (has_column(base, base_sel))
                base = base.Filter([](bool pass) { return pass; }, {base_sel});
            else
                base = base.Filter(base_sel);
        }

        auto mask_ext = el.mask_for_ext();
        auto mask_mc_like = el.mask_for_mc_like();

        // All non-data samples.
        ROOT::RDF::RNode node_all =
            filter_by_sample_mask(base, mask_mc_like, "sample_id");

        // MC only, i.e. non-EXT non-data.
        ROOT::RDF::RNode node_mc =
            filter_not_sample_mask(node_all, mask_ext, "sample_id");

        // Signal is MC only.
        ROOT::RDF::RNode node_sig = node_mc.Filter(signal_sel);

        // Total background = MC background (+ EXT if requested).
        ROOT::RDF::RNode node_bkg =
            include_ext_in_background
                ? node_all.Filter("!(" + signal_sel + ")")
                : node_mc.Filter("!(" + signal_sel + ")");

        const ULong64_t n_sig_raw = *node_sig.Count();
        const ULong64_t n_bkg_raw = *node_bkg.Count();

        if (n_sig_raw == 0 && n_bkg_raw == 0)
        {
            std::cerr << "[plot_nominal_weight_signal_bkg_unstacked] no selected events after base_sel='"
                      << base_sel << "'.\n";
            return 1;
        }

        // Auto-range if no valid [xmin, xmax] was supplied.
        if (!(xmax > xmin))
        {
            double lo = std::numeric_limits<double>::infinity();
            double hi = -std::numeric_limits<double>::infinity();

            if (n_sig_raw > 0)
            {
                const double smin = *node_sig.Min<double>("__plot_x__");
                const double smax = *node_sig.Max<double>("__plot_x__");
                lo = std::min(lo, smin);
                hi = std::max(hi, smax);
            }

            if (n_bkg_raw > 0)
            {
                const double bmin = *node_bkg.Min<double>("__plot_x__");
                const double bmax = *node_bkg.Max<double>("__plot_x__");
                lo = std::min(lo, bmin);
                hi = std::max(hi, bmax);
            }

            if (!std::isfinite(lo) || !std::isfinite(hi))
            {
                std::cerr << "[plot_nominal_weight_signal_bkg_unstacked] could not determine a valid x-range.\n";
                return 1;
            }

            if (hi <= lo)
            {
                const double eps = (std::abs(lo) > 0.0) ? 0.05 * std::abs(lo) : 1.0;
                lo -= eps;
                hi += eps;
            }

            const double pad = 0.05 * (hi - lo);
            xmin = lo - pad;
            xmax = hi + pad;
        }

        ROOT::RDF::TH1DModel h_sig_model("h_nominal_weight_sig_raw", "", nbins, xmin, xmax);
        ROOT::RDF::TH1DModel h_bkg_model("h_nominal_weight_bkg_raw", "", nbins, xmin, xmax);

        // Unweighted fills: x-axis is weight_col itself.
        auto h_sig = node_sig.Histo1D(h_sig_model, "__plot_x__");
        auto h_bkg = node_bkg.Histo1D(h_bkg_model, "__plot_x__");

        TH1D &hs = *h_sig;
        TH1D &hb = *h_bkg;

        set_hist_style(hs, kRed + 1, 20);
        set_hist_style(hb, kBlue + 1, 24);

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        TCanvas c("c_nominal_weight_signal_bkg_unstacked",
                  "Nominal weight: signal vs background",
                  1100,
                  800);
        c.SetLeftMargin(0.12);
        c.SetRightMargin(0.04);
        c.SetBottomMargin(0.12);
        if (use_logy)
            c.SetLogy();

        const double ymax = std::max(hs.GetMaximum(), hb.GetMaximum());
        const double frame_max = use_logy
                                     ? 10.0 * std::max(1.0, ymax)
                                     : 1.25 * std::max(1.0, ymax);

        const std::string frame_title =
            ";" + weight_col + ";Raw event count";

        TH1D h_frame("h_frame_nominal_weight", frame_title.c_str(), 100, xmin, xmax);
        h_frame.SetMinimum(use_logy ? 0.5 : 0.0);
        h_frame.SetMaximum(frame_max);
        h_frame.Draw("AXIS");

        hb.Draw("HIST SAME");
        hs.Draw("HIST SAME");

        TLegend leg(0.56, 0.72, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);

        std::ostringstream sig_lbl;
        sig_lbl << "signal (MC only), N = " << n_sig_raw;

        std::ostringstream bkg_lbl;
        if (include_ext_in_background)
            bkg_lbl << "total background (MC bkg + EXT), N = " << n_bkg_raw;
        else
            bkg_lbl << "background (MC only), N = " << n_bkg_raw;

        leg.AddEntry(&hb, bkg_lbl.str().c_str(), "l");
        leg.AddEntry(&hs, sig_lbl.str().c_str(), "l");
        leg.Draw();

        c.RedrawAxis();

        const auto out = plot_output_file(output_stem).string();
        c.SaveAs(out.c_str());

        std::cout << "[plot_nominal_weight_signal_bkg_unstacked] base_sel='"
                  << base_sel << "' signal_sel='" << signal_sel << "'\n";
        std::cout << "  signal raw count = " << n_sig_raw << "\n";
        std::cout << "  bkg raw count    = " << n_bkg_raw << "\n";
        std::cout << "  x-range          = [" << xmin << ", " << xmax << "]\n";
        std::cout << "[plot_nominal_weight_signal_bkg_unstacked] saved: "
                  << out << "\n";

        return 0;
    });
}
