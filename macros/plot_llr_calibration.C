// plot_llr_calibration.C
//
// Diagnostic plots to test whether inf_scores[0] behaves like a log-likelihood ratio
// between signal and background.
//
// The macro builds score histograms for truth signal (is_signal) and truth background
// (!is_signal) after a base selection, then constructs two binned tests:
//
//  (A) "LLR in score space" (priors removed):
//        y(s) = log p(s | signal) - log p(s | background)
//      where p(s|class) is estimated from the per-class normalised score histogram
//      (optionally weighted). If inf_scores[0] is the true LLR, then y(s) = s.
//
//  (B) "Posterior log-odds" under the sample priors:
//        y_post(s) = log P(signal | s-bin) - log P(background | s-bin)
//                  = log N_sig(bin) - log N_bkg(bin)
//      If inf_scores[0] is the true LLR, then y_post(s) = s + log(pi_sig/pi_bkg),
//      with pi_sig/pi_bkg estimated from the total (weighted) yields.
//
// Uncertainties:
//   Points use per-bin statistical uncertainties from the (weighted) histograms,
//   propagated to the log-ratio:  sigma^2 = (sigma_sig/n_sig)^2 + (sigma_bkg/n_bkg)^2.
//   (Normalisation uncertainties are neglected; for large totals this is subdominant.)
//
// Run with:
//   ./heron macro plot_llr_calibration.C
//   ./heron macro plot_llr_calibration.C \
//     'plot_llr_calibration("./scratch/out/event_list_myana.root")'
//
// Example with custom binning/range:
//   ./heron macro plot_llr_calibration.C \
//     'plot_llr_calibration("./scratch/out/event_list_myana.root", "sel_muon", "is_signal", "w_nominal", 60, -15, 15, true, "first_score_llr_check")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
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
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

bool implicit_mt_enabled()
{
    const char *env = std::getenv("HERON_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

double safe_log(double x)
{
    // Caller should ensure x>0; this is a last-resort guard.
    if (!(x > 0.0) || !std::isfinite(x))
        return std::numeric_limits<double>::quiet_NaN();
    return std::log(x);
}

struct Point
{
    double x = 0.0;
    double ex = 0.0;
    double y = 0.0;
    double ey = 0.0;
};

struct CalibCurve
{
    std::vector<double> x;
    std::vector<double> y;
};

CalibCurve isotonic_curve(const std::vector<Point> &pts)
{
    struct Block
    {
        int first = 0;
        int last = 0;
        double w = 1.0;
        double mean = 0.0;
    };

    std::vector<Block> blocks;
    blocks.reserve(pts.size());

    for (int i = 0; i < static_cast<int>(pts.size()); ++i)
    {
        double w = 1.0;
        if (pts[i].ey > 0.0 && std::isfinite(pts[i].ey))
            w = 1.0 / (pts[i].ey * pts[i].ey);

        blocks.push_back(Block{i, i, w, pts[i].y});

        while (blocks.size() >= 2)
        {
            auto &a = blocks[blocks.size() - 2];
            auto &b = blocks[blocks.size() - 1];

            if (a.mean <= b.mean)
                break;

            const double wtot = a.w + b.w;
            const double mtot = (a.w * a.mean + b.w * b.mean) / wtot;
            a = Block{a.first, b.last, wtot, mtot};
            blocks.pop_back();
        }
    }

    CalibCurve out;
    out.x.reserve(pts.size());
    out.y.reserve(pts.size());

    for (const auto &blk : blocks)
    {
        for (int i = blk.first; i <= blk.last; ++i)
        {
            out.x.push_back(pts[i].x);
            out.y.push_back(blk.mean);
        }
    }

    return out;
}

double eval_piecewise_linear(const std::vector<double> &x,
                             const std::vector<double> &y,
                             double s)
{
    if (x.empty())
        return s;

    if (s <= x.front())
        return y.front();
    if (s >= x.back())
        return y.back();

    auto it = std::upper_bound(x.begin(), x.end(), s);
    const std::size_t i1 = static_cast<std::size_t>(it - x.begin());
    const std::size_t i0 = i1 - 1;

    const double dx = x[i1] - x[i0];
    if (!(dx > 0.0))
        return y[i0];

    const double t = (s - x[i0]) / dx;
    return y[i0] + t * (y[i1] - y[i0]);
}

std::vector<Point> build_ratio_points(const TH1D &hs,
                                      const TH1D &hb,
                                      double tot_sig,
                                      double tot_bkg,
                                      bool posterior)
{
    std::vector<Point> pts;
    const int nbins = hs.GetNbinsX();
    pts.reserve(static_cast<std::size_t>(nbins));

    for (int b = 1; b <= nbins; ++b)
    {
        const double ns = hs.GetBinContent(b);
        const double nb = hb.GetBinContent(b);
        const double ens = hs.GetBinError(b);
        const double enb = hb.GetBinError(b);

        if (!(ns > 0.0) || !(nb > 0.0))
            continue;

        const double width = hs.GetBinWidth(b);
        if (!(width > 0.0))
            continue;

        const double x = hs.GetBinCenter(b);
        const double ex = 0.5 * width;

        double y = 0.0;
        if (posterior)
            y = safe_log(ns) - safe_log(nb);
        else
            y = safe_log(ns / (tot_sig * width)) - safe_log(nb / (tot_bkg * width));

        double ey = 0.0;
        if (ens > 0.0 && enb > 0.0)
        {
            const double rs = ens / ns;
            const double rb = enb / nb;
            ey = std::sqrt(rs * rs + rb * rb);
        }

        pts.push_back(Point{x, ex, y, ey});
    }

    return pts;
}

} // namespace

int plotFirstInferenceScoreLLRCalibration(const std::string &event_list_path = "",
                                         const std::string &base_sel = "sel_muon",
                                         const std::string &signal_sel = "is_signal",
                                         const std::string &mc_weight = "w_nominal",
                                         int nbins = 60,
                                         double xmin = -15.0,
                                         double xmax = 15.0,
                                         bool use_weights = true,
                                         const std::string &output_stem = "first_inference_score_llr_calibration")
{
    return heron::macro::run_with_guard("plotFirstInferenceScoreLLRCalibration", [&]() -> int {
        if (nbins < 1)
            nbins = 1;
        if (xmax < xmin)
            std::swap(xmin, xmax);
        if (xmax == xmin)
            xmax = xmin + 1.0;

        if (implicit_mt_enabled())
            ROOT::EnableImplicitMT();

        TH1::SetDefaultSumw2();

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plotFirstInferenceScoreLLRCalibration] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plotFirstInferenceScoreLLRCalibration] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode rdf = SelectionService::decorate(el.rdf())
                               .Define("inf_score_0",
                                       [](const ROOT::RVec<float> &scores) {
                                           return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                       },
                                       {"inf_scores"});

        if (!has_column(rdf, signal_sel))
        {
            std::cerr << "[plotFirstInferenceScoreLLRCalibration] missing signal selection column/expression: '"
                      << signal_sel << "'.\n";
            std::cerr << "  Expected a boolean column like 'is_signal' in the event list.\n";
            return 1;
        }

        auto mask_mc_like = el.mask_for_mc_like();
        ROOT::RDF::RNode node = filter_by_sample_mask(rdf, mask_mc_like, "sample_id")
                                    .Define("__w__", mc_weight);

        if (!base_sel.empty())
        {
            if (has_column(rdf, base_sel))
                node = node.Filter([](bool pass) { return pass; }, {base_sel});
            else
                node = node.Filter(base_sel);
        }

        node = node.Filter(
            [xmin, xmax](double s) {
                return std::isfinite(s) && (s >= xmin) && (s < xmax);
            },
            {"inf_score_0"});

        ROOT::RDF::RNode node_calib = node.Filter("rdfentry_ % 2 == 0");
        ROOT::RDF::RNode node_valid = node.Filter("rdfentry_ % 2 == 1");

        ROOT::RDF::RNode node_sig = node_calib.Filter(signal_sel);
        ROOT::RDF::RNode node_bkg = node_calib.Filter("!(" + signal_sel + ")");

        ROOT::RDF::TH1DModel hmodel_sig("h_sig", "", nbins, xmin, xmax);
        ROOT::RDF::TH1DModel hmodel_bkg("h_bkg", "", nbins, xmin, xmax);

        ROOT::RDF::RResultPtr<TH1D> h_sig = use_weights
                                                ? node_sig.Histo1D(hmodel_sig, "inf_score_0", "__w__")
                                                : node_sig.Histo1D(hmodel_sig, "inf_score_0");

        ROOT::RDF::RResultPtr<TH1D> h_bkg = use_weights
                                                ? node_bkg.Histo1D(hmodel_bkg, "inf_score_0", "__w__")
                                                : node_bkg.Histo1D(hmodel_bkg, "inf_score_0");

        const TH1D &hs = *h_sig;
        const TH1D &hb = *h_bkg;

        const double tot_sig = hs.Integral(1, nbins);
        const double tot_bkg = hb.Integral(1, nbins);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[plotFirstInferenceScoreLLRCalibration] totals in range [" << xmin << ", " << xmax << "):\n";
        std::cout << "  signal total=" << tot_sig << "\n";
        std::cout << "  bkg    total=" << tot_bkg << "\n";

        if (!(tot_sig > 0.0) || !(tot_bkg > 0.0))
        {
            std::cerr << "[plotFirstInferenceScoreLLRCalibration] need positive signal and background totals in the plotted range.\n";
            return 1;
        }

        const double prior_offset = std::log(tot_sig / tot_bkg);
        std::cout << "[plotFirstInferenceScoreLLRCalibration] log(prior odds) ~ log(tot_sig/tot_bkg) = "
                  << prior_offset << "\n";

        std::vector<Point> pts_llr = build_ratio_points(hs, hb, tot_sig, tot_bkg, false);
        std::vector<Point> pts_post = build_ratio_points(hs, hb, tot_sig, tot_bkg, true);

        if (pts_llr.size() < 2)
        {
            std::cerr << "[plotFirstInferenceScoreLLRCalibration] not enough populated bins with both signal and background to build the LLR check graph.\n";
            std::cerr << "  Try fewer/wider bins or a wider [xmin,xmax] range.\n";
            return 1;
        }

        const auto curve = isotonic_curve(pts_llr);
        const auto xcal = curve.x;
        const auto ycal = curve.y;

        ROOT::RDF::RNode node_valid_cal = node_valid.Define(
            "inf_score_0_cal",
            [xcal, ycal](double s) { return eval_piecewise_linear(xcal, ycal, s); },
            {"inf_score_0"});

        ROOT::RDF::RNode node_valid_sig = node_valid_cal.Filter(signal_sel);
        ROOT::RDF::RNode node_valid_bkg = node_valid_cal.Filter("!(" + signal_sel + ")");

        ROOT::RDF::TH1DModel hmodel_sig_cal("h_sig_cal", "", nbins, xmin, xmax);
        ROOT::RDF::TH1DModel hmodel_bkg_cal("h_bkg_cal", "", nbins, xmin, xmax);

        ROOT::RDF::RResultPtr<TH1D> h_sig_cal = use_weights
                                                    ? node_valid_sig.Histo1D(hmodel_sig_cal, "inf_score_0_cal", "__w__")
                                                    : node_valid_sig.Histo1D(hmodel_sig_cal, "inf_score_0_cal");

        ROOT::RDF::RResultPtr<TH1D> h_bkg_cal = use_weights
                                                    ? node_valid_bkg.Histo1D(hmodel_bkg_cal, "inf_score_0_cal", "__w__")
                                                    : node_valid_bkg.Histo1D(hmodel_bkg_cal, "inf_score_0_cal");

        const TH1D &hs_cal = *h_sig_cal;
        const TH1D &hb_cal = *h_bkg_cal;
        const double tot_sig_cal = hs_cal.Integral(1, nbins);
        const double tot_bkg_cal = hb_cal.Integral(1, nbins);

        std::vector<Point> pts_llr_cal;
        if (tot_sig_cal > 0.0 && tot_bkg_cal > 0.0)
            pts_llr_cal = build_ratio_points(hs_cal, hb_cal, tot_sig_cal, tot_bkg_cal, false);

        double y_llr_min = std::numeric_limits<double>::infinity();
        double y_llr_max = -std::numeric_limits<double>::infinity();
        double y_post_min = std::numeric_limits<double>::infinity();
        double y_post_max = -std::numeric_limits<double>::infinity();

        for (const auto &p : pts_llr)
        {
            y_llr_min = std::min(y_llr_min, p.y);
            y_llr_max = std::max(y_llr_max, p.y);
        }
        for (const auto &p : pts_post)
        {
            y_post_min = std::min(y_post_min, p.y);
            y_post_max = std::max(y_post_max, p.y);
        }

        auto make_graph = [](const std::vector<Point> &pts, const char *name) {
            std::vector<double> x, y, ex, ey;
            x.reserve(pts.size());
            y.reserve(pts.size());
            ex.reserve(pts.size());
            ey.reserve(pts.size());
            for (const auto &p : pts)
            {
                x.push_back(p.x);
                y.push_back(p.y);
                ex.push_back(p.ex);
                ey.push_back(p.ey);
            }
            auto g = std::make_unique<TGraphErrors>(static_cast<int>(pts.size()), x.data(), y.data(), ex.data(), ey.data());
            g->SetName(name);
            g->SetLineWidth(2);
            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.0);
            return g;
        };

        auto g_llr = make_graph(pts_llr, "g_llr_score_space");
        auto g_post = make_graph(pts_post, "g_post_logodds");
        std::unique_ptr<TGraphErrors> g_llr_cal;
        if (pts_llr_cal.size() >= 2)
            g_llr_cal = make_graph(pts_llr_cal, "g_llr_score_space_calibrated");

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        {
            TCanvas c("c_score_llr_check", "inf_scores[0] LLR check", 1100, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            if (!std::isfinite(y_llr_min) || !std::isfinite(y_llr_max))
            {
                y_llr_min = xmin;
                y_llr_max = xmax;
            }
            const double ypad = 0.10 * std::max(1.0, y_llr_max - y_llr_min);
            const double y_min = y_llr_min - ypad;
            const double y_max = y_llr_max + ypad;

            TH1D hframe("hframe_llr",
                        ";Inference score [0];log p(score|signal) - log p(score|background)",
                        100,
                        xmin,
                        xmax);
            hframe.SetMinimum(y_min);
            hframe.SetMaximum(y_max);
            hframe.Draw("AXIS");

            g_llr->Draw("P SAME");
            g_llr->Draw("L SAME");

            TF1 f_id("f_id", "x", xmin, xmax);
            f_id.SetLineStyle(2);
            f_id.SetLineWidth(2);
            f_id.Draw("SAME");

            TF1 f_fit("f_fit", "[0] + [1]*x", xmin, xmax);
            f_fit.SetLineStyle(1);
            f_fit.SetLineWidth(2);
            g_llr->Fit(&f_fit, "Q");

            std::cout << "[plotFirstInferenceScoreLLRCalibration] score-space LLR fit: y = p0 + p1*x\n";
            std::cout << "  p0 = " << f_fit.GetParameter(0) << " +- " << f_fit.GetParError(0) << "\n";
            std::cout << "  p1 = " << f_fit.GetParameter(1) << " +- " << f_fit.GetParError(1) << "\n";

            TLegend leg(0.13, 0.73, 0.62, 0.89);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_llr.get(), "binned log p(s|S) - log p(s|B)", "lp");
            leg.AddEntry(&f_id, "y = x (ideal LLR)", "l");
            leg.AddEntry(&f_fit, "linear fit", "l");
            leg.Draw();

            c.RedrawAxis();

            const auto out = plot_output_file(output_stem + "_score_space_llr").string();
            c.SaveAs(out.c_str());
            std::cout << "[plotFirstInferenceScoreLLRCalibration] saved: " << out << "\n";
        }

        if (g_llr_cal)
        {
            TCanvas c("c_score_llr_check_calibrated", "Calibrated inf_scores[0] LLR check", 1100, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            double ycal_min = std::numeric_limits<double>::infinity();
            double ycal_max = -std::numeric_limits<double>::infinity();
            for (const auto &p : pts_llr_cal)
            {
                ycal_min = std::min(ycal_min, p.y);
                ycal_max = std::max(ycal_max, p.y);
            }
            if (!std::isfinite(ycal_min) || !std::isfinite(ycal_max))
            {
                ycal_min = xmin;
                ycal_max = xmax;
            }

            const double ypad = 0.10 * std::max(1.0, ycal_max - ycal_min);
            TH1D hframe("hframe_llr_cal",
                        ";Calibrated inference score [0];log p(score|signal) - log p(score|background)",
                        100,
                        xmin,
                        xmax);
            hframe.SetMinimum(ycal_min - ypad);
            hframe.SetMaximum(ycal_max + ypad);
            hframe.Draw("AXIS");

            g_llr_cal->SetMarkerStyle(22);
            g_llr_cal->Draw("P SAME");
            g_llr_cal->Draw("L SAME");

            TF1 f_id("f_id_cal", "x", xmin, xmax);
            f_id.SetLineStyle(2);
            f_id.SetLineWidth(2);
            f_id.Draw("SAME");

            TF1 f_fit("f_fit_cal", "[0] + [1]*x", xmin, xmax);
            f_fit.SetLineStyle(1);
            f_fit.SetLineWidth(2);
            g_llr_cal->Fit(&f_fit, "Q");

            std::cout << "[plotFirstInferenceScoreLLRCalibration] calibrated score-space LLR fit (validation): y = p0 + p1*x\n";
            std::cout << "  p0 = " << f_fit.GetParameter(0) << " +- " << f_fit.GetParError(0) << "\n";
            std::cout << "  p1 = " << f_fit.GetParameter(1) << " +- " << f_fit.GetParError(1) << "\n";

            TLegend leg(0.13, 0.73, 0.69, 0.89);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_llr_cal.get(), "validation: binned log p(s_{cal}|S) - log p(s_{cal}|B)", "lp");
            leg.AddEntry(&f_id, "y = x (ideal closure)", "l");
            leg.AddEntry(&f_fit, "linear fit", "l");
            leg.Draw();

            c.RedrawAxis();

            const auto out = plot_output_file(output_stem + "_score_space_llr_calibrated").string();
            c.SaveAs(out.c_str());
            std::cout << "[plotFirstInferenceScoreLLRCalibration] saved: " << out << "\n";
        }
        else
        {
            std::cout << "[plotFirstInferenceScoreLLRCalibration] skipped calibrated closure plot: not enough validation bins with both classes.\n";
        }

        {
            TCanvas c("c_score_post_logodds", "inf_scores[0] posterior log-odds", 1100, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            if (!std::isfinite(y_post_min) || !std::isfinite(y_post_max))
            {
                y_post_min = xmin;
                y_post_max = xmax;
            }
            const double ypad = 0.10 * std::max(1.0, y_post_max - y_post_min);
            const double y_min = y_post_min - ypad;
            const double y_max = y_post_max + ypad;

            TH1D hframe("hframe_post",
                        ";Inference score [0];log P(signal|bin) - log P(background|bin)",
                        100,
                        xmin,
                        xmax);
            hframe.SetMinimum(y_min);
            hframe.SetMaximum(y_max);
            hframe.Draw("AXIS");

            g_post->SetMarkerStyle(21);
            g_post->Draw("P SAME");
            g_post->Draw("L SAME");

            TF1 f_expect("f_expect", "x + [0]", xmin, xmax);
            f_expect.SetParameter(0, prior_offset);
            f_expect.SetLineStyle(2);
            f_expect.SetLineWidth(2);
            f_expect.Draw("SAME");

            TF1 f_fit("f_fit_post", "[0] + [1]*x", xmin, xmax);
            f_fit.SetLineStyle(1);
            f_fit.SetLineWidth(2);
            g_post->Fit(&f_fit, "Q");

            std::cout << "[plotFirstInferenceScoreLLRCalibration] posterior log-odds fit: y = p0 + p1*x\n";
            std::cout << "  p0 = " << f_fit.GetParameter(0) << " +- " << f_fit.GetParError(0) << "\n";
            std::cout << "  p1 = " << f_fit.GetParameter(1) << " +- " << f_fit.GetParError(1) << "\n";

            std::ostringstream prior_ss;
            prior_ss << "y = x + log(totS/totB) = x + " << std::setprecision(3) << prior_offset;

            TLegend leg(0.13, 0.73, 0.72, 0.89);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_post.get(), "binned posterior log-odds (raw bin counts)", "lp");
            leg.AddEntry(&f_expect, prior_ss.str().c_str(), "l");
            leg.AddEntry(&f_fit, "linear fit", "l");
            leg.Draw();

            c.RedrawAxis();

            const auto out = plot_output_file(output_stem + "_posterior_logodds").string();
            c.SaveAs(out.c_str());
            std::cout << "[plotFirstInferenceScoreLLRCalibration] saved: " << out << "\n";
        }

        std::cout << "[plotFirstInferenceScoreLLRCalibration] done\n";
        return 0;
    });
}

int plot_llr_calibration(const std::string &event_list_path = "",
                         const std::string &base_sel = "sel_muon",
                         const std::string &signal_sel = "is_signal",
                         const std::string &mc_weight = "w_nominal",
                         int nbins = 60,
                         double xmin = -15.0,
                         double xmax = 15.0,
                         bool use_weights = true,
                         const std::string &output_stem = "plot_llr_calibration")
{
    return plotFirstInferenceScoreLLRCalibration(event_list_path, base_sel, signal_sel, mc_weight, nbins,
                                                  xmin, xmax, use_weights, output_stem);
}
