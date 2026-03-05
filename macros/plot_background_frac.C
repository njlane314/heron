// plot/macro/plot_background_frac.C
//
// Companion to plot_background_rejection.C.
//
// Scan a threshold on inf_scores[0] and, for each background analysis channel,
// plot the background passing fraction
//
//   epsilon_pass(thr) = N_pass(thr) / N_total
//
// on a log-y axis.
//
// Statistical treatment:
//   - N_total and N_pass are RAW event counts after the requested base
//     selection.
//   - Pointwise uncertainty bars are exact Clopper-Pearson intervals on the
//     pass fraction.
//   - Weighted MC/EXT yields are still reported in the printed table at the
//     chosen report threshold, but the plotted intervals are based on raw counts.
//
// Important detail for a log-y plot:
//   - Points with N_pass = 0 cannot be drawn at y = 0 on a logarithmic axis.
//   - Once a channel first reaches zero passing raw events, this macro draws a
//     dashed horizontal line at the chosen one-sided upper limit for the pass
//     fraction (default: 95% C.L.), starting from that threshold. An open marker
//     is placed at the transition point.
//   - These dashed segments are upper limits, not central-value measurements.
//
// Suggested usage:
//   ./heron macro plot_background_frac.C
//   ./heron macro plot_background_frac.C \
//     'plot_background_frac("./scratch/out/event_list_myana.root")'
//   ./heron macro plot_background_frac.C \
//     'plot_background_frac("./scratch/out/event_list_myana.root", "sel_muon", "w_nominal", 81, -15., 15., 4.0)'

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

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

struct ChannelDraw
{
    ChannelStyle style;
    ULong64_t n_total_raw = 0;
    std::unique_ptr<TGraphAsymmErrors> g_pass;
    std::unique_ptr<TGraph> g_ul_marker;
    std::unique_ptr<TLine> l_ul;
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

        out.push_back(ChannelStyle{key, channel_colour(key), marker});

        marker += 1;
        if (marker > 34)
            marker = 20;
    }

    out.push_back(ChannelStyle{
        static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown),
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

void set_pass_graph_style(TGraphAsymmErrors &g, const ChannelStyle &style)
{
    g.SetLineColor(style.color);
    g.SetMarkerColor(style.color);
    g.SetLineWidth(2);
    g.SetMarkerStyle(style.marker);
    g.SetMarkerSize(1.0);
}

void set_ul_marker_style(TGraph &g, const ChannelStyle &style)
{
    g.SetLineColor(style.color);
    g.SetMarkerColor(style.color);
    g.SetLineWidth(2);
    g.SetMarkerStyle(24);
    g.SetMarkerSize(1.1);
}

void set_ul_line_style(TLine &l, const ChannelStyle &style)
{
    l.SetLineColor(style.color);
    l.SetLineStyle(2);
    l.SetLineWidth(2);
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

double positive_axis_floor(double value)
{
    if (!(value > 0.0) || !std::isfinite(value))
        return 1.0e-6;

    double floor = std::pow(10.0, std::floor(std::log10(value)));
    if (!(floor > 0.0) || !std::isfinite(floor))
        floor = value;

    if (floor >= value)
        floor *= 0.5;

    return std::max(1.0e-12, floor);
}

} // namespace

int plot_background_frac(const std::string &event_list_path = "",
                         const std::string &base_sel = "sel_muon",
                         const std::string &mc_weight = "w_nominal",
                         int n_thresholds = 81,
                         double raw_threshold_min = -15.0,
                         double raw_threshold_max = 15.0,
                         double report_threshold = 4.0,
                         double cl = 0.682689492137086,
                         double zero_ul_cl = 0.954499736103642,
                         const std::string &output_stem = "first_inference_score_background_passfrac_log")
{
    return heron::macro::run_with_guard("plot_background_frac", [&]() -> int {
        ROOT::EnableImplicitMT();
        TH1::SetDefaultSumw2();

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_background_frac] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_background_frac] input is not an event-list root file: "
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
            std::cerr << "[plot_background_frac] missing column 'analysis_channels'.\n"
                      << "  The event list needs to include the derived truth-level analysis channel column.\n";
            return 1;
        }

        auto mask_bkg = el.mask_for_mc_like();

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
            ROOT::RDF::TH1DModel hmodel(("h_raw_" + suff).c_str(), "", n_thresholds, edges.data());

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
                node_ch.Histo1D(hmodel, "inf_score_0"),
                node_report.Count(),
                node_report.Sum<double>("__w__"),
                node_report.Sum<double>("__w2__")});
        }

        std::vector<ChannelDraw> drawings;
        std::vector<ChannelReport> reports;
        drawings.reserve(booked.size());
        reports.reserve(booked.size());

        double min_positive_shown = std::numeric_limits<double>::infinity();

        for (auto &book : booked)
        {
            const ULong64_t n_total_raw = *book.n_total_raw;
            const double w_total = *book.w_total;
            const double w2_total = *book.w2_total;

            if (n_total_raw == 0)
                continue;

            const TH1D &h_raw = *book.h_raw;

            std::vector<double> x_pass;
            std::vector<double> y_pass;
            std::vector<double> exl_pass;
            std::vector<double> exh_pass;
            std::vector<double> eyl_pass;
            std::vector<double> eyh_pass;

            bool have_zero_region = false;
            double first_zero_threshold = raw_threshold_max;
            double zero_ul = pass_fraction_hi_cp(n_total_raw, 0, zero_ul_cl);
            if (!(zero_ul > 0.0) || !std::isfinite(zero_ul))
                zero_ul = 1.0 / static_cast<double>(std::max<ULong64_t>(1, n_total_raw + 1ULL));

            for (int i = 0; i < n_thresholds; ++i)
            {
                const int first_bin = i + 1;
                const ULong64_t n_pass_raw = static_cast<ULong64_t>(std::llround(cumulative_hist_sum(h_raw, first_bin)));
                const double thr = thresholds[static_cast<std::size_t>(i)];

                if (n_pass_raw > 0)
                {
                    const double pass = static_cast<double>(n_pass_raw) / static_cast<double>(n_total_raw);
                    const double pass_lo = pass_fraction_lo_cp(n_total_raw, n_pass_raw, cl);
                    const double pass_hi = pass_fraction_hi_cp(n_total_raw, n_pass_raw, cl);

                    x_pass.push_back(thr);
                    y_pass.push_back(pass);
                    exl_pass.push_back(0.0);
                    exh_pass.push_back(0.0);
                    eyl_pass.push_back(pass - pass_lo);
                    eyh_pass.push_back(pass_hi - pass);

                    min_positive_shown = std::min(min_positive_shown, pass_lo);
                }
                else if (!have_zero_region)
                {
                    have_zero_region = true;
                    first_zero_threshold = thr;
                    min_positive_shown = std::min(min_positive_shown, zero_ul);
                }
            }

            ChannelDraw draw;
            draw.style = book.style;
            draw.n_total_raw = n_total_raw;

            if (!x_pass.empty())
            {
                draw.g_pass = std::make_unique<TGraphAsymmErrors>(static_cast<int>(x_pass.size()),
                                                                  x_pass.data(),
                                                                  y_pass.data(),
                                                                  exl_pass.data(),
                                                                  exh_pass.data(),
                                                                  eyl_pass.data(),
                                                                  eyh_pass.data());
                set_pass_graph_style(*draw.g_pass, draw.style);
            }

            if (have_zero_region)
            {
                double x_ul[1] = {first_zero_threshold};
                double y_ul[1] = {zero_ul};
                draw.g_ul_marker = std::make_unique<TGraph>(1, x_ul, y_ul);
                set_ul_marker_style(*draw.g_ul_marker, draw.style);

                draw.l_ul = std::make_unique<TLine>(first_zero_threshold, zero_ul,
                                                    raw_threshold_max, zero_ul);
                set_ul_line_style(*draw.l_ul, draw.style);
            }

            drawings.push_back(std::move(draw));

            reports.push_back(ChannelReport{
                book.style,
                n_total_raw,
                *book.n_report_raw,
                w_total,
                w2_total,
                *book.w_report,
                *book.w2_report});
        }

        if (drawings.empty())
        {
            std::cerr << "[plot_background_frac] no non-empty background channels after base selection='"
                      << base_sel << "'.\n";
            return 1;
        }

        if (!std::isfinite(min_positive_shown))
            min_positive_shown = 1.0e-6;

        const double y_min = positive_axis_floor(min_positive_shown);
        const double y_max = 1.3;

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        TCanvas c("c_plot_background_frac", "Background pass fraction by analysis channel", 1150, 800);
        c.SetLeftMargin(0.11);
        c.SetRightMargin(0.04);
        c.SetBottomMargin(0.12);
        c.SetLogy();

        TH1D h_frame("h_frame_bkgpasslog",
                     ";Inference score [0] threshold;background pass fraction #varepsilon_{pass}",
                     100,
                     raw_threshold_min,
                     raw_threshold_max);
        h_frame.SetMinimum(y_min);
        h_frame.SetMaximum(y_max);
        h_frame.Draw("AXIS");

        for (auto &d : drawings)
        {
            if (d.g_pass)
            {
                d.g_pass->Draw("PZ SAME");
                d.g_pass->Draw("L SAME");
            }
            if (d.l_ul)
                d.l_ul->Draw("SAME");
            if (d.g_ul_marker)
                d.g_ul_marker->Draw("P SAME");
        }

        TLine report_line(report_threshold, y_min, report_threshold, y_max);
        report_line.SetLineStyle(2);
        report_line.SetLineWidth(2);
        report_line.Draw();

        double dummy_x[1] = {0.0};
        double dummy_y[1] = {1.0};
        TGraph dummy_ul(1, dummy_x, dummy_y);
        dummy_ul.SetMarkerStyle(24);
        dummy_ul.SetMarkerSize(1.1);
        dummy_ul.SetLineStyle(2);
        dummy_ul.SetLineWidth(2);

        TLegend leg(0.13, 0.14, 0.96, 0.42);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetNColumns(3);
        for (auto &d : drawings)
        {
            std::ostringstream entry;
            entry << channel_label(d.style.id) << " (N=" << d.n_total_raw << ")";
            if (d.g_pass)
                leg.AddEntry(d.g_pass.get(), entry.str().c_str(), "lp");
            else if (d.g_ul_marker)
                leg.AddEntry(d.g_ul_marker.get(), entry.str().c_str(), "lp");
            else if (d.l_ul)
                leg.AddEntry(d.l_ul.get(), entry.str().c_str(), "l");
        }
        leg.AddEntry(&dummy_ul, "open marker + dashed line: 95% UL after first 0-pass point", "lp");
        std::ostringstream report_label;
        report_label << "report threshold = " << std::setprecision(3) << report_threshold;
        leg.AddEntry(&report_line, report_label.str().c_str(), "l");
        leg.Draw();

        c.RedrawAxis();

        const auto out = plot_output_file(output_stem).string();
        c.SaveAs(out.c_str());

        const double cl95 = 0.954499736103642;
        std::cout << "\n[plot_background_frac] report at threshold = "
                  << report_threshold << "\n";
        std::cout << "  base_sel = '" << base_sel << "'\n";
        std::cout << "  raw-count intervals are exact Clopper-Pearson on the pass fraction\n";
        std::cout << "  weighted uncertainties are sumw2 statistical errors on the surviving yield\n";
        std::cout << "  in the log-y plot, 0-pass regions are shown as 95% upper-limit segments\n\n";

        std::cout << std::left
                  << std::setw(26) << "channel"
                  << std::setw(14) << "raw pass/tot"
                  << std::setw(24) << "pass frac [68% C.L.]"
                  << std::setw(18) << "pass frac 95% UL"
                  << std::setw(14) << "rejection"
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

            std::ostringstream pass68_ss;
            pass68_ss << std::fixed << std::setprecision(4)
                      << pass << " [" << pass_lo68 << ", " << pass_hi68 << "]";

            std::ostringstream pass95_ss;
            pass95_ss << std::fixed << std::setprecision(4) << pass_hi95;

            std::ostringstream rej_ss;
            rej_ss << std::fixed << std::setprecision(4) << rejection;

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
                      << std::setw(24) << pass68_ss.str()
                      << std::setw(18) << pass95_ss.str()
                      << std::setw(14) << rej_ss.str()
                      << std::setw(22) << wpass_ss.str()
                      << std::setw(14) << wtot_ss.str()
                      << std::setw(12) << nefftot_ss.str()
                      << neffpass_ss.str()
                      << "\n";
        }

        std::cout << "\n[plot_background_frac] saved plot: " << out << "\n";
        std::cout << "[plot_background_frac] done\n";

        return 0;
    });
}
