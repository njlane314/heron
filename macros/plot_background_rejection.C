// plot/macro/plot_background_rejection.C
//
// Scan a threshold on inf_scores[0] and, for each background analysis channel,
// plot the background rejection
//
//   rejection(thr) = 1 - N_pass(thr) / N_total
//
// where N_total and N_pass are RAW event counts after the requested base
// selection. The uncertainty bars are exact Clopper-Pearson binomial intervals
// on the passing fraction, mapped onto the rejection axis.
//
// Additional plots produced by this macro:
//   1) Pass-fraction (raw) on a log-y axis:
//        eps_pass(thr)      = N_pass(thr) / N_total
//        eps_pass^95%UL(thr)= CP-upper(N_total, N_pass(thr); CL=0.95)
//      The 95% UL is overlaid as a dashed curve (same colour).
//
//   2) Effective weighted tail statistics vs threshold:
//        N_eff,pass(thr) = (sum w_pass)^2 / sum (w_pass^2)
//      This indicates whether a small weighted tail is supported by meaningful
//      statistics in the high-score region.
//
// Why raw counts for the bars?
//   - Exact binomial intervals are well-defined for unweighted counts.
//   - For weighted MC/EXT events, the ratio of weighted sums is not strictly a
//     binomial quantity.
//   - To address the "is there really zero background in the signal region?"
//     question, the macro also prints a per-channel report at a chosen threshold
//     with:
//       * raw pass / total counts
//       * 68% and 95% upper limits on the pass fraction
//       * weighted surviving yield and its sumw2 statistical uncertainty
//       * effective weighted statistics N_eff = (sum w)^2 / sum w^2
//
// Suggested usage:
//   ./heron macro plot_background_rejection.C
//   ./heron macro plot_background_rejection.C \
//     'plot_background_rejection("./scratch/out/event_list_myana.root")'
//   ./heron macro plot_background_rejection.C \
//     'plot_background_rejection("./scratch/out/event_list_myana.root", "sel_muon", "w_nominal", 81, -15., 15., 4.0)'
//
// Notes:
//   - The uncertainty bars are pointwise. Adjacent threshold points are highly
//     correlated because the scan is cumulative.
//   - If you want the visually more informative version near the signal region,
//     keep this rejection plot and additionally inspect the printed table at the
//     chosen report_threshold.

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
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
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

#include "AnalysisChannels.hh"
#include "EventListIO.hh"
#include "PlotChannels.hh"
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

struct ChannelStyle
{
    int id = 0;
    const char *label = "";
    int color = kBlack;
    int marker = 20;
};

struct ChannelActions
{
    ChannelStyle style;
    ROOT::RDF::RResultPtr<ULong64_t> n_total_raw;
    ROOT::RDF::RResultPtr<double> w_total;
    ROOT::RDF::RResultPtr<double> w2_total;
    ROOT::RDF::RResultPtr<TH1D> h_raw;
    ROOT::RDF::RResultPtr<TH1D> h_w;
    ROOT::RDF::RResultPtr<TH1D> h_w2;
    ROOT::RDF::RResultPtr<ULong64_t> n_report_raw;
    ROOT::RDF::RResultPtr<double> w_report;
    ROOT::RDF::RResultPtr<double> w2_report;
};

struct ChannelReport
{
    ChannelStyle style;
    ULong64_t n_total_raw = 0;
    ULong64_t n_pass_raw = 0;
    double w_total = 0.0;
    double w2_total = 0.0;
    double w_pass = 0.0;
    double w2_pass = 0.0;
};

std::string channel_label(int id)
{
    if (id == static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown))
        return "Unknown";
    return Channels::label(id);
}

int channel_colour(int id)
{
    if (id == static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown))
        return kBlack;
    return Channels::colour(id);
}

std::vector<ChannelStyle> default_background_channels()
{
    std::vector<ChannelStyle> out;
    int marker = 20;

    for (int key : Channels::mc_keys())
    {
        if (key == static_cast<int>(AnalysisChannels::AnalysisChannel::SignalLambda))
            continue;
        if (key == static_cast<int>(AnalysisChannels::AnalysisChannel::DataInclusive))
            continue;

        out.push_back(ChannelStyle{
            key,
            nullptr,
            channel_colour(key),
            marker});

        marker += 1;
        if (marker > 34)
            marker = 20;
    }

    out.push_back(ChannelStyle{
        static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown),
        nullptr,
        channel_colour(static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown)),
        25});

    return out;
}

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

std::string channel_suffix(int id)
{
    std::ostringstream ss;
    ss << "ch" << id;
    return ss.str();
}

void set_graph_style(TGraphAsymmErrors &g, const ChannelStyle &style)
{
    g.SetLineColor(style.color);
    g.SetMarkerColor(style.color);
    g.SetLineWidth(2);
    g.SetMarkerStyle(style.marker);
    g.SetMarkerSize(1.0);
}

void set_graph_style(TGraph &g, const ChannelStyle &style)
{
    g.SetLineColor(style.color);
    g.SetMarkerColor(style.color);
    g.SetLineWidth(2);
    g.SetMarkerStyle(style.marker);
    g.SetMarkerSize(1.0);
}

double cumulative_hist_sum(const TH1 &h, int first_bin)
{
    const int last_bin = h.GetNbinsX() + 1; // include overflow
    double out = 0.0;
    for (int b = first_bin; b <= last_bin; ++b)
        out += h.GetBinContent(b);
    return out;
}

double pass_fraction_lo_cp(ULong64_t n_total, ULong64_t n_pass, double cl)
{
    return TEfficiency::ClopperPearson(static_cast<unsigned int>(n_total),
                                       static_cast<unsigned int>(n_pass),
                                       cl,
                                       false);
}

double pass_fraction_hi_cp(ULong64_t n_total, ULong64_t n_pass, double cl)
{
    return TEfficiency::ClopperPearson(static_cast<unsigned int>(n_total),
                                       static_cast<unsigned int>(n_pass),
                                       cl,
                                       true);
}

double effective_entries(double sumw, double sumw2)
{
    if (sumw2 <= 0.0)
        return 0.0;
    return (sumw * sumw) / sumw2;
}

} // namespace

int plot_background_rejection(const std::string &event_list_path = "",
                                               const std::string &base_sel = "sel_muon",
                                               const std::string &mc_weight = "w_nominal",
                                               int n_thresholds = 81,
                                               double raw_threshold_min = -15.0,
                                               double raw_threshold_max = 15.0,
                                               double report_threshold = 4.0,
                                               double cl = 0.682689492137086,
                                               const std::string &output_stem = "plot_background_rejection")
{
    return heron::macro::run_with_guard("plot_background_rejection", [&]() -> int {
        ROOT::EnableImplicitMT();
        TH1::SetDefaultSumw2();

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_background_rejection] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_background_rejection] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        if (n_thresholds < 2)
            n_thresholds = 2;

        if (raw_threshold_max < raw_threshold_min)
            std::swap(raw_threshold_min, raw_threshold_max);
        if (raw_threshold_max == raw_threshold_min)
            raw_threshold_max = raw_threshold_min + 1.0;

        const double step = (raw_threshold_max - raw_threshold_min) / static_cast<double>(n_thresholds - 1);

        std::vector<double> thresholds(static_cast<std::size_t>(n_thresholds));
        std::vector<double> edges(static_cast<std::size_t>(n_thresholds + 1));
        for (int i = 0; i < n_thresholds; ++i)
        {
            const double thr = raw_threshold_min + static_cast<double>(i) * step;
            thresholds[static_cast<std::size_t>(i)] = thr;
            edges[static_cast<std::size_t>(i)] = thr;
        }
        edges.back() = raw_threshold_max + step;

        EventListIO el(input_path);

        ROOT::RDF::RNode rdf = SelectionService::decorate(el.rdf())
                                   .Define("inf_score_0",
                                           [](const ROOT::RVec<float> &scores) {
                                               return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                           },
                                           {"inf_scores"});

        if (!has_column(rdf, "analysis_channels"))
        {
            std::cerr << "[plot_background_rejection] missing column 'analysis_channels'.\n"
                      << "  The event list needs to include the derived truth-level analysis channel column.\n";
            return 1;
        }

        auto mask_bkg = el.mask_for_mc_like(); // non-data: MC + EXT

        ROOT::RDF::RNode node_bkg = filter_by_sample_mask(rdf, mask_bkg, "sample_id")
                                        .Define("__w__", mc_weight)
                                        .Define("__w2__", "__w__ * __w__");

        if (!base_sel.empty())
        {
            if (has_column(rdf, base_sel))
            {
                node_bkg = node_bkg.Filter([](bool pass) { return pass; }, {base_sel});
            }
            else
            {
                node_bkg = node_bkg.Filter(base_sel);
            }
        }

        const auto channels = default_background_channels();
        std::vector<ChannelActions> booked;
        booked.reserve(channels.size());

        for (const auto &style : channels)
        {
            ROOT::RDF::RNode node_ch = node_bkg.Filter(
                [cid = style.id](int channel) {
                    return channel == cid;
                },
                {"analysis_channels"});

            const std::string suff = channel_suffix(style.id);
            ROOT::RDF::TH1DModel hmodel_raw(("h_raw_" + suff).c_str(), "", n_thresholds, edges.data());
            ROOT::RDF::TH1DModel hmodel_w(("h_w_" + suff).c_str(), "", n_thresholds, edges.data());
            ROOT::RDF::TH1DModel hmodel_w2(("h_w2_" + suff).c_str(), "", n_thresholds, edges.data());

            ROOT::RDF::RNode node_report = node_ch.Filter(
                [report_threshold](double score) {
                    return score >= report_threshold;
                },
                {"inf_score_0"});

            booked.push_back(ChannelActions{
                style,
                node_ch.Count(),
                node_ch.Sum<double>("__w__"),
                node_ch.Sum<double>("__w2__"),
                node_ch.Histo1D(hmodel_raw, "inf_score_0"),
                node_ch.Histo1D(hmodel_w, "inf_score_0", "__w__"),
                node_ch.Histo1D(hmodel_w2, "inf_score_0", "__w2__"),
                node_report.Count(),
                node_report.Sum<double>("__w__"),
                node_report.Sum<double>("__w2__")});
        }

        // --- Per-channel graphs ---
        std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs_rej;
        std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs_pass;
        std::vector<std::unique_ptr<TGraph>> graphs_pass_ul95;
        std::vector<std::unique_ptr<TGraph>> graphs_neff_pass;
        std::vector<ChannelReport> reports;
        graphs_rej.reserve(booked.size());
        graphs_pass.reserve(booked.size());
        graphs_pass_ul95.reserve(booked.size());
        graphs_neff_pass.reserve(booked.size());
        reports.reserve(booked.size());

        // Global ranges for the added plots.
        double min_pass_like = std::numeric_limits<double>::infinity();
        double max_neff_pass = 0.0;

        // 95% confidence level for the CP upper-limit overlay.
        const double cl95 = 0.95;

        for (auto &book : booked)
        {
            const ULong64_t n_total_raw = *book.n_total_raw;
            const double w_total = *book.w_total;
            const double w2_total = *book.w2_total;

            if (n_total_raw == 0)
                continue;

            const TH1D &h_raw = *book.h_raw;
            const TH1D &h_w = *book.h_w;
            const TH1D &h_w2 = *book.h_w2;

            std::vector<double> x(static_cast<std::size_t>(n_thresholds));
            std::vector<double> y(static_cast<std::size_t>(n_thresholds));
            std::vector<double> exl(static_cast<std::size_t>(n_thresholds), 0.0);
            std::vector<double> exh(static_cast<std::size_t>(n_thresholds), 0.0);
            std::vector<double> eyl(static_cast<std::size_t>(n_thresholds), 0.0);
            std::vector<double> eyh(static_cast<std::size_t>(n_thresholds), 0.0);

            // Pass-fraction graph (skip y=0 points for log-y plotting).
            std::vector<double> x_pass;
            std::vector<double> y_pass;
            std::vector<double> exl_pass;
            std::vector<double> exh_pass;
            std::vector<double> eyl_pass;
            std::vector<double> eyh_pass;
            x_pass.reserve(static_cast<std::size_t>(n_thresholds));
            y_pass.reserve(static_cast<std::size_t>(n_thresholds));
            exl_pass.reserve(static_cast<std::size_t>(n_thresholds));
            exh_pass.reserve(static_cast<std::size_t>(n_thresholds));
            eyl_pass.reserve(static_cast<std::size_t>(n_thresholds));
            eyh_pass.reserve(static_cast<std::size_t>(n_thresholds));

            // 95% UL curve (always positive for n_total_raw > 0).
            std::vector<double> x_ul95(static_cast<std::size_t>(n_thresholds));
            std::vector<double> y_ul95(static_cast<std::size_t>(n_thresholds));

            // N_eff,pass curve (skip 0 to allow optional log-y plotting).
            std::vector<double> x_neff;
            std::vector<double> y_neff;
            x_neff.reserve(static_cast<std::size_t>(n_thresholds));
            y_neff.reserve(static_cast<std::size_t>(n_thresholds));

            for (int i = 0; i < n_thresholds; ++i)
            {
                const int first_bin = i + 1;
                const ULong64_t n_pass_raw = static_cast<ULong64_t>(std::llround(cumulative_hist_sum(h_raw, first_bin)));
                const double pass = static_cast<double>(n_pass_raw) / static_cast<double>(n_total_raw);
                const double pass_lo = pass_fraction_lo_cp(n_total_raw, n_pass_raw, cl);
                const double pass_hi = pass_fraction_hi_cp(n_total_raw, n_pass_raw, cl);
                const double rejection = 1.0 - pass;

                x[static_cast<std::size_t>(i)] = thresholds[static_cast<std::size_t>(i)];
                y[static_cast<std::size_t>(i)] = rejection;

                // Convert [pass_lo, pass_hi] into the corresponding rejection interval
                // [1 - pass_hi, 1 - pass_lo]. The asymmetric errors therefore swap.
                eyl[static_cast<std::size_t>(i)] = pass_hi - pass;
                eyh[static_cast<std::size_t>(i)] = pass - pass_lo;

                // --- Added plot 1: pass fraction (log-y) with 95% UL overlay ---
                const double pass_ul = pass_fraction_hi_cp(n_total_raw, n_pass_raw, cl95);
                x_ul95[static_cast<std::size_t>(i)] = thresholds[static_cast<std::size_t>(i)];
                y_ul95[static_cast<std::size_t>(i)] = pass_ul;
                if (pass_ul > 0.0)
                    min_pass_like = std::min(min_pass_like, pass_ul);
                if (pass > 0.0)
                {
                    x_pass.push_back(thresholds[static_cast<std::size_t>(i)]);
                    y_pass.push_back(pass);
                    exl_pass.push_back(0.0);
                    exh_pass.push_back(0.0);
                    eyl_pass.push_back(pass - pass_lo);
                    eyh_pass.push_back(pass_hi - pass);
                    min_pass_like = std::min(min_pass_like, pass);
                }

                // --- Added plot 2: N_eff,pass(thr) from cumulative weighted sums ---
                const double sumw_pass = cumulative_hist_sum(h_w, first_bin);
                const double sumw2_pass = cumulative_hist_sum(h_w2, first_bin);
                const double neff_pass = effective_entries(sumw_pass, sumw2_pass);
                if (neff_pass > 0.0)
                {
                    x_neff.push_back(thresholds[static_cast<std::size_t>(i)]);
                    y_neff.push_back(neff_pass);
                    max_neff_pass = std::max(max_neff_pass, neff_pass);
                }
            }

            // Rejection graph.
            auto g_rej = std::make_unique<TGraphAsymmErrors>(n_thresholds,
                                                             x.data(),
                                                             y.data(),
                                                             exl.data(),
                                                             exh.data(),
                                                             eyl.data(),
                                                             eyh.data());
            set_graph_style(*g_rej, book.style);
            graphs_rej.push_back(std::move(g_rej));

            // Pass-fraction graph (raw), with CP interval error bars at `cl`.
            auto g_pass = std::make_unique<TGraphAsymmErrors>(static_cast<int>(x_pass.size()),
                                                              x_pass.data(),
                                                              y_pass.data(),
                                                              exl_pass.data(),
                                                              exh_pass.data(),
                                                              eyl_pass.data(),
                                                              eyh_pass.data());
            set_graph_style(*g_pass, book.style);
            graphs_pass.push_back(std::move(g_pass));

            // Pass-fraction 95% UL curve (dashed, same colour).
            auto g_ul = std::make_unique<TGraph>(n_thresholds, x_ul95.data(), y_ul95.data());
            set_graph_style(*g_ul, book.style);
            g_ul->SetLineStyle(2);
            g_ul->SetMarkerSize(0.0);
            graphs_pass_ul95.push_back(std::move(g_ul));

            // N_eff,pass curve.
            auto g_neff = std::make_unique<TGraph>(static_cast<int>(x_neff.size()), x_neff.data(), y_neff.data());
            set_graph_style(*g_neff, book.style);
            graphs_neff_pass.push_back(std::move(g_neff));

            reports.push_back(ChannelReport{
                book.style,
                n_total_raw,
                *book.n_report_raw,
                w_total,
                w2_total,
                *book.w_report,
                *book.w2_report});
        }

        if (graphs_rej.empty())
        {
            std::cerr << "[plot_background_rejection] no non-empty background channels after base selection='"
                      << base_sel << "'.\n";
            return 1;
        }

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);
        TCanvas c("c_first_inf_score_bkgrej", "Background rejection by analysis channel", 1150, 800);
        c.SetLeftMargin(0.11);
        c.SetRightMargin(0.04);
        c.SetBottomMargin(0.12);

        TH1D h_frame("h_frame_bkgrej",
                     ";Inference score [0] threshold;background rejection = 1 - #varepsilon_{pass}",
                     100,
                     raw_threshold_min,
                     raw_threshold_max);
        h_frame.SetMinimum(0.0);
        h_frame.SetMaximum(1.05);
        h_frame.Draw("AXIS");

        for (auto &gptr : graphs_rej)
        {
            gptr->Draw("PZ SAME");
            gptr->Draw("L SAME");
        }

        TLine report_line(report_threshold, 0.0, report_threshold, 1.05);
        report_line.SetLineStyle(2);
        report_line.SetLineWidth(2);
        report_line.Draw();

        TLegend leg(0.13, 0.14, 0.96, 0.40);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetNColumns(3);
        for (std::size_t i = 0; i < graphs_rej.size(); ++i)
        {
            const auto &r = reports[i];
            std::ostringstream entry;
            entry << channel_label(r.style.id) << " (N=" << r.n_total_raw << ")";
            leg.AddEntry(graphs_rej[i].get(), entry.str().c_str(), "lp");
        }
        std::ostringstream report_label;
        report_label << "report threshold = " << std::setprecision(3) << report_threshold;
        leg.AddEntry(&report_line, report_label.str().c_str(), "l");
        leg.Draw();

        c.RedrawAxis();

        const auto out = plot_output_file(output_stem).string();
        c.SaveAs(out.c_str());

        // --------------------------------------------------------------------
        // Added plot 1: pass fraction on log-y with 95% UL dashed curve.
        // --------------------------------------------------------------------
        {
            plotter.set_global_style();
            TCanvas c_pass("c_first_inf_score_bkgpass", "Background pass fraction by analysis channel", 1150, 800);
            c_pass.SetLeftMargin(0.11);
            c_pass.SetRightMargin(0.04);
            c_pass.SetBottomMargin(0.12);
            c_pass.SetLogy();

            if (!std::isfinite(min_pass_like) || min_pass_like <= 0.0)
                min_pass_like = 1.0e-8;

            // Pad minimum: slightly below the smallest UL/positive point.
            const double y_min_pass = std::max(1.0e-8, 0.5 * min_pass_like);

            TH1D h_frame_pass("h_frame_bkgpass",
                              ";Inference score [0] threshold;background pass fraction #varepsilon_{pass}",
                              100,
                              raw_threshold_min,
                              raw_threshold_max);
            h_frame_pass.SetMinimum(y_min_pass);
            h_frame_pass.SetMaximum(1.05);
            h_frame_pass.Draw("AXIS");

            // Draw central values first, then the UL overlay (dashed).
            for (std::size_t i = 0; i < graphs_pass.size(); ++i)
            {
                graphs_pass[i]->Draw("PZ SAME");
                graphs_pass[i]->Draw("L SAME");
                graphs_pass_ul95[i]->Draw("L SAME");
            }

            TLine report_line_pass(report_threshold, y_min_pass, report_threshold, 1.05);
            report_line_pass.SetLineStyle(2);
            report_line_pass.SetLineWidth(2);
            report_line_pass.Draw();

            // Style legend: keep channel colours + a short style key.
            TLine solid_key(raw_threshold_min, 1.0, raw_threshold_min + 0.1, 1.0);
            solid_key.SetLineStyle(1);
            solid_key.SetLineWidth(2);
            TLine dashed_key(raw_threshold_min, 1.0, raw_threshold_min + 0.1, 1.0);
            dashed_key.SetLineStyle(2);
            dashed_key.SetLineWidth(2);

            TLegend leg_pass(0.13, 0.14, 0.96, 0.43);
            leg_pass.SetBorderSize(0);
            leg_pass.SetFillStyle(0);
            leg_pass.SetNColumns(3);
            leg_pass.AddEntry(&solid_key, "#hat{#varepsilon}_{pass} (raw)", "l");
            leg_pass.AddEntry(&dashed_key, "95% UL (Clopper-Pearson)", "l");
            leg_pass.AddEntry(&report_line_pass, report_label.str().c_str(), "l");

            for (std::size_t i = 0; i < graphs_pass.size(); ++i)
            {
                const auto &r = reports[i];
                std::ostringstream entry;
                entry << channel_label(r.style.id) << " (N=" << r.n_total_raw << ")";
                leg_pass.AddEntry(graphs_pass[i].get(), entry.str().c_str(), "lp");
            }
            leg_pass.Draw();

            c_pass.RedrawAxis();
            const auto out_pass = plot_output_file(output_stem + "_pass_log").string();
            c_pass.SaveAs(out_pass.c_str());
            std::cout << "[plot_background_rejection] saved pass-fraction plot: " << out_pass << "\n";
        }

        // --------------------------------------------------------------------
        // Added plot 2: N_eff,pass(threshold) from weighted passing tail.
        // --------------------------------------------------------------------
        {
            plotter.set_global_style();
            TCanvas c_neff("c_first_inf_score_bkgneff", "Effective passing statistics by analysis channel", 1150, 800);
            c_neff.SetLeftMargin(0.11);
            c_neff.SetRightMargin(0.04);
            c_neff.SetBottomMargin(0.12);
            c_neff.SetLogy();

            const double y_min_neff = 0.5;
            const double y_max_neff = (max_neff_pass > 0.0) ? (1.2 * max_neff_pass) : 10.0;

            TH1D h_frame_neff("h_frame_bkgneff",
                              ";Inference score [0] threshold;N_{eff,pass} = (#Sigma w)^{2}/(#Sigma w^{2})",
                              100,
                              raw_threshold_min,
                              raw_threshold_max);
            h_frame_neff.SetMinimum(y_min_neff);
            h_frame_neff.SetMaximum(y_max_neff);
            h_frame_neff.Draw("AXIS");

            for (std::size_t i = 0; i < graphs_neff_pass.size(); ++i)
            {
                graphs_neff_pass[i]->Draw("LP SAME");
            }

            TLine report_line_neff(report_threshold, y_min_neff, report_threshold, y_max_neff);
            report_line_neff.SetLineStyle(2);
            report_line_neff.SetLineWidth(2);
            report_line_neff.Draw();

            TLegend leg_neff(0.13, 0.14, 0.96, 0.40);
            leg_neff.SetBorderSize(0);
            leg_neff.SetFillStyle(0);
            leg_neff.SetNColumns(3);
            for (std::size_t i = 0; i < graphs_neff_pass.size(); ++i)
            {
                const auto &r = reports[i];
                std::ostringstream entry;
                entry << channel_label(r.style.id) << " (N=" << r.n_total_raw << ")";
                leg_neff.AddEntry(graphs_neff_pass[i].get(), entry.str().c_str(), "lp");
            }
            leg_neff.AddEntry(&report_line_neff, report_label.str().c_str(), "l");
            leg_neff.Draw();

            c_neff.RedrawAxis();
            const auto out_neff = plot_output_file(output_stem + "_neff_pass").string();
            c_neff.SaveAs(out_neff.c_str());
            std::cout << "[plot_background_rejection] saved N_eff,pass plot: " << out_neff << "\n";
        }

        std::cout << "\n[plot_background_rejection] report at threshold = "
                  << report_threshold << "\n";
        std::cout << "  base_sel = '" << base_sel << "'\n";
        std::cout << "  raw-count intervals are exact Clopper-Pearson on the pass fraction\n";
        std::cout << "  weighted uncertainties are sumw2 statistical errors on the surviving yield\n\n";

        std::cout << std::left
                  << std::setw(26) << "channel"
                  << std::setw(14) << "raw pass/tot"
                  << std::setw(14) << "rej"
                  << std::setw(30) << "pass frac [68% C.L.]"
                  << std::setw(18) << "pass frac 95% UL"
                  << std::setw(22) << "weighted pass +- stat"
                  << std::setw(14) << "w total"
                  << std::setw(12) << "N_eff tot"
                  << "N_eff pass"
                  << "\n";

        for (const auto &r : reports)
        {
            const double pass = (r.n_total_raw > 0)
                                    ? static_cast<double>(r.n_pass_raw) / static_cast<double>(r.n_total_raw)
                                    : 0.0;
            const double pass_lo68 = pass_fraction_lo_cp(r.n_total_raw, r.n_pass_raw, cl);
            const double pass_hi68 = pass_fraction_hi_cp(r.n_total_raw, r.n_pass_raw, cl);
            const double pass_hi95 = pass_fraction_hi_cp(r.n_total_raw, r.n_pass_raw, cl95);
            const double rejection = 1.0 - pass;
            const double w_pass_err = std::sqrt(std::max(0.0, r.w2_pass));
            const double n_eff_tot = effective_entries(r.w_total, r.w2_total);
            const double n_eff_pass = effective_entries(r.w_pass, r.w2_pass);

            std::ostringstream raw_ratio;
            raw_ratio << r.n_pass_raw << "/" << r.n_total_raw;

            std::ostringstream rej_ss;
            rej_ss << std::fixed << std::setprecision(4) << rejection;

            std::ostringstream pass68_ss;
            pass68_ss << std::fixed << std::setprecision(4)
                      << pass << " [" << pass_lo68 << ", " << pass_hi68 << "]";

            std::ostringstream pass95_ss;
            pass95_ss << std::fixed << std::setprecision(4) << pass_hi95;

            std::ostringstream wpass_ss;
            wpass_ss << std::fixed << std::setprecision(3)
                     << r.w_pass << " +- " << w_pass_err;

            std::ostringstream wtot_ss;
            wtot_ss << std::fixed << std::setprecision(3) << r.w_total;

            std::ostringstream nefftot_ss;
            nefftot_ss << std::fixed << std::setprecision(1) << n_eff_tot;

            std::ostringstream neffpass_ss;
            neffpass_ss << std::fixed << std::setprecision(1) << n_eff_pass;

            std::cout << std::left
                      << std::setw(26) << channel_label(r.style.id)
                      << std::setw(14) << raw_ratio.str()
                      << std::setw(14) << rej_ss.str()
                      << std::setw(30) << pass68_ss.str()
                      << std::setw(18) << pass95_ss.str()
                      << std::setw(22) << wpass_ss.str()
                      << std::setw(14) << wtot_ss.str()
                      << std::setw(12) << nefftot_ss.str()
                      << neffpass_ss.str()
                      << "\n";
        }

        std::cout << "\n[plot_background_rejection] saved plot: " << out << "\n";
        std::cout << "[plot_background_rejection] done\n";

        return 0;
    });
}

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif
