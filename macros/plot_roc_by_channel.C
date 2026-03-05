// plot/macro/plot_roc_by_channel.C
//
// Scan a threshold on inf_scores[0] and build ROC-like curves per analysis
// channel.
//
// For each threshold t:
//   x(t) = signal efficiency          = N_S,pass(t) / N_S,total
//   y_c(t)= background rejection (c)  = 1 - N_c,pass(t) / N_c,total
//
// where all N are RAW event counts after the requested base selection.
//
// Uncertainties:
//   - exact Clopper-Pearson intervals at confidence level `cl` are used for
//     both the signal efficiency and the background passing fraction.
//   - the background interval is mapped onto the rejection axis.
//
// Notes:
//   - This scan is cumulative in threshold, so adjacent points are highly
//     correlated.
//   - The CP intervals are strictly defined for unweighted counts; this is a
//     feature here (you see whether you have enough raw tail statistics).
//   - The macro also draws an "All backgrounds" curve (everything with
//     !(signal_sel)).
//
// Run with:
//   ./heron macro plot_roc_by_channel.C
//   ./heron macro plot_roc_by_channel.C \
//     'plot_roc_by_channel("./scratch/out/event_list_myana.root")'
//   ./heron macro plot_roc_by_channel.C \
//     'plot_roc_by_channel("./scratch/out/event_list_myana.root", "sel_muon", "is_signal", 81, -15, 15)'

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
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

struct ChannelStyle
{
    int id = 0;
    std::string label;
    int color = kBlack;
    int marker = 20;
};

struct ChannelBooked
{
    ChannelStyle style;
    ROOT::RDF::RResultPtr<ULong64_t> n_total_raw;
    ROOT::RDF::RResultPtr<TH1D> h_raw;
};

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

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

        out.push_back(ChannelStyle{key, channel_label(key), channel_colour(key), marker});
        marker += 1;
        if (marker > 34)
            marker = 20;
    }

    out.push_back(ChannelStyle{static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown),
                               channel_label(static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown)),
                               channel_colour(static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown)),
                               25});

    return out;
}

void set_graph_style(TGraphAsymmErrors &g, const ChannelStyle &style, int lw = 2)
{
    g.SetLineColor(style.color);
    g.SetMarkerColor(style.color);
    g.SetLineWidth(lw);
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

double cp_lo(ULong64_t n_total, ULong64_t n_pass, double cl)
{
    return TEfficiency::ClopperPearson(static_cast<unsigned int>(n_total),
                                       static_cast<unsigned int>(n_pass),
                                       cl,
                                       false);
}

double cp_hi(ULong64_t n_total, ULong64_t n_pass, double cl)
{
    return TEfficiency::ClopperPearson(static_cast<unsigned int>(n_total),
                                       static_cast<unsigned int>(n_pass),
                                       cl,
                                       true);
}

double auc_trapezoid(const std::vector<double> &x, const std::vector<double> &y)
{
    if (x.size() < 2 || y.size() != x.size())
        return 0.0;

    // Sort by x ascending.
    std::vector<std::size_t> idx(x.size());
    for (std::size_t i = 0; i < idx.size(); ++i)
        idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });

    double area = 0.0;
    for (std::size_t i = 0; i + 1 < idx.size(); ++i)
    {
        const std::size_t ia = idx[i];
        const std::size_t ib = idx[i + 1];
        const double dx = x[ib] - x[ia];
        area += 0.5 * (y[ia] + y[ib]) * dx;
    }
    return area;
}

} // namespace

int plot_roc_by_channel(const std::string &event_list_path = "",
                        const std::string &base_sel = "sel_muon",
                        const std::string &signal_sel = "is_signal",
                        int n_thresholds = 81,
                        double raw_threshold_min = -15.0,
                        double raw_threshold_max = 15.0,
                        double cl = 0.682689492137086,
                        const std::string &output_stem = "plot_roc_by_channel")
{
    return heron::macro::run_with_guard("plot_roc_by_channel", [&]() -> int {
        ROOT::EnableImplicitMT();
        TH1::SetDefaultSumw2(false);

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_roc_by_channel] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_roc_by_channel] input is not an event-list root file: "
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
            std::cerr << "[plot_roc_by_channel] missing column 'analysis_channels'.\n";
            return 1;
        }

        auto mask_mc_like = el.mask_for_mc_like();
        ROOT::RDF::RNode node = filter_by_sample_mask(rdf, mask_mc_like, "sample_id");

        if (!base_sel.empty())
        {
            if (has_column(rdf, base_sel))
                node = node.Filter([](bool pass) { return pass; }, {base_sel});
            else
                node = node.Filter(base_sel);
        }

        // Signal booking.
        ROOT::RDF::RNode node_sig = node.Filter(signal_sel);
        ROOT::RDF::TH1DModel hmodel_sig("h_sig_raw", "", n_thresholds, edges.data());
        auto n_sig_total = node_sig.Count();
        auto h_sig_raw = node_sig.Histo1D(hmodel_sig, "inf_score_0");

        // All-background booking.
        ROOT::RDF::RNode node_bkg_all = node.Filter("!(" + signal_sel + ")");
        ROOT::RDF::TH1DModel hmodel_bkg_all("h_bkg_all_raw", "", n_thresholds, edges.data());
        auto n_bkg_all_total = node_bkg_all.Count();
        auto h_bkg_all_raw = node_bkg_all.Histo1D(hmodel_bkg_all, "inf_score_0");

        // Per-channel background booking.
        std::vector<ChannelBooked> booked;
        const auto channels = default_background_channels();
        booked.reserve(channels.size());

        for (const auto &style : channels)
        {
            ROOT::RDF::RNode node_ch = node.Filter(
                [cid = style.id](int ch) { return ch == cid; },
                {"analysis_channels"});

            ROOT::RDF::TH1DModel hmodel(("h_bkg_" + std::to_string(style.id)).c_str(), "", n_thresholds, edges.data());

            booked.push_back(ChannelBooked{style, node_ch.Count(), node_ch.Histo1D(hmodel, "inf_score_0")});
        }

        const ULong64_t Nsig = *n_sig_total;
        if (Nsig == 0)
        {
            std::cerr << "[plot_roc_by_channel] no signal events after base selection='" << base_sel
                      << "' and signal_sel='" << signal_sel << "'.\n";
            return 1;
        }

        const TH1D &hS = *h_sig_raw;

        // Precompute signal efficiency and CP errors for each threshold.
        std::vector<double> x_eff(static_cast<std::size_t>(n_thresholds));
        std::vector<double> x_el(static_cast<std::size_t>(n_thresholds));
        std::vector<double> x_eh(static_cast<std::size_t>(n_thresholds));

        for (int i = 0; i < n_thresholds; ++i)
        {
            const int first_bin = i + 1;
            const ULong64_t n_pass = static_cast<ULong64_t>(std::llround(cumulative_hist_sum(hS, first_bin)));
            const double eff = static_cast<double>(n_pass) / static_cast<double>(Nsig);
            const double lo = cp_lo(Nsig, n_pass, cl);
            const double hi = cp_hi(Nsig, n_pass, cl);

            x_eff[static_cast<std::size_t>(i)] = eff;
            x_el[static_cast<std::size_t>(i)] = eff - lo;
            x_eh[static_cast<std::size_t>(i)] = hi - eff;
        }

        // Make a ROC graph from a background histogram.
        auto make_roc = [&](const TH1D &hB, ULong64_t Nbkg, const ChannelStyle &style, int lw) {
            std::vector<double> y_rej(static_cast<std::size_t>(n_thresholds));
            std::vector<double> y_el(static_cast<std::size_t>(n_thresholds));
            std::vector<double> y_eh(static_cast<std::size_t>(n_thresholds));
            std::vector<double> exl(static_cast<std::size_t>(n_thresholds), 0.0);
            std::vector<double> exh(static_cast<std::size_t>(n_thresholds), 0.0);

            for (int i = 0; i < n_thresholds; ++i)
            {
                const int first_bin = i + 1;
                const ULong64_t n_pass = static_cast<ULong64_t>(std::llround(cumulative_hist_sum(hB, first_bin)));
                const double pass = static_cast<double>(n_pass) / static_cast<double>(Nbkg);
                const double lo = cp_lo(Nbkg, n_pass, cl);
                const double hi = cp_hi(Nbkg, n_pass, cl);

                const double rej = 1.0 - pass;
                y_rej[static_cast<std::size_t>(i)] = rej;

                // Map pass interval [lo, hi] to rejection interval [1-hi, 1-lo]
                y_el[static_cast<std::size_t>(i)] = hi - pass; // lower error on rejection
                y_eh[static_cast<std::size_t>(i)] = pass - lo; // upper error on rejection
            }

            auto g = std::make_unique<TGraphAsymmErrors>(n_thresholds,
                                                         x_eff.data(),
                                                         y_rej.data(),
                                                         x_el.data(),
                                                         x_eh.data(),
                                                         y_el.data(),
                                                         y_eh.data());
            set_graph_style(*g, style, lw);
            return g;
        };

        std::vector<std::unique_ptr<TGraphAsymmErrors>> rocs;
        std::vector<ChannelStyle> roc_styles;

        // All backgrounds first (for axes).
        const ULong64_t NbkgAll = n_bkg_all_total.GetValue();
        if (NbkgAll > 0)
        {
            ChannelStyle all_style;
            all_style.id = -1;
            all_style.label = "All backgrounds";
            all_style.color = kBlack;
            all_style.marker = 20;
            rocs.push_back(make_roc(h_bkg_all_raw.GetValue(), NbkgAll, all_style, 3));
            roc_styles.push_back(all_style);
        }

        // Per-channel.
        for (const auto &b : booked)
        {
            const ULong64_t Nbkg = b.n_total_raw.GetValue();
            if (Nbkg == 0)
                continue;
            // Skip the signal channel explicitly (defensive).
            if (b.style.id == static_cast<int>(AnalysisChannels::AnalysisChannel::SignalLambda))
                continue;

            rocs.push_back(make_roc(b.h_raw.GetValue(), Nbkg, b.style, 2));
            roc_styles.push_back(b.style);
        }

        if (rocs.empty())
        {
            std::cerr << "[plot_roc_by_channel] no backgrounds after base selection='" << base_sel << "'.\n";
            return 1;
        }

        // Print AUCs.
        std::cout << "\n[plot_roc_by_channel] AUC( signal eff vs background rejection )\n";
        for (std::size_t i = 0; i < rocs.size(); ++i)
        {
            const int n = rocs[i]->GetN();
            std::vector<double> x(n), y(n);
            for (int j = 0; j < n; ++j)
            {
                double xx = 0.0, yy = 0.0;
                rocs[i]->GetPoint(j, xx, yy);
                x[static_cast<std::size_t>(j)] = xx;
                y[static_cast<std::size_t>(j)] = yy;
            }
            const double auc = auc_trapezoid(x, y);
            std::cout << "  " << std::setw(24) << roc_styles[i].label << " : "
                      << std::fixed << std::setprecision(4) << auc << "\n";
        }

        // Plot.
        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        TCanvas c("c_roc_by_channel", "ROC by analysis channel", 1100, 850);
        c.SetLeftMargin(0.11);
        c.SetRightMargin(0.04);
        c.SetBottomMargin(0.12);

        // Axis frame.
        TH1D h_frame("h_frame_roc", ";signal efficiency #varepsilon_{S};background rejection 1 - #varepsilon_{B}", 100, 0.0, 1.0);
        h_frame.SetMinimum(0.0);
        h_frame.SetMaximum(1.05);
        h_frame.Draw("AXIS");

        for (auto &g : rocs)
        {
            g->Draw("PZ SAME");
            g->Draw("L SAME");
        }

        // Reference lines.
        TLine diag(0.0, 0.0, 1.0, 1.0);
        diag.SetLineStyle(2);
        diag.SetLineWidth(1);
        diag.SetLineColor(16);
        diag.Draw();

        TLegend leg(0.12, 0.14, 0.96, 0.40);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetNColumns(3);

        for (std::size_t i = 0; i < rocs.size(); ++i)
        {
            std::ostringstream ss;
            ss << roc_styles[i].label;
            // Add totals for non-all backgrounds when available.
            if (roc_styles[i].id >= 0)
            {
                // Find the matching booked entry to get N.
                ULong64_t N = 0;
                for (const auto &b : booked)
                {
                    if (b.style.id == roc_styles[i].id)
                    {
                        N = b.n_total_raw.GetValue();
                        break;
                    }
                }
                if (N > 0)
                    ss << " (N=" << N << ")";
            }
            else
            {
                ss << " (N=" << NbkgAll << ")";
            }
            leg.AddEntry(rocs[i].get(), ss.str().c_str(), "lp");
        }
        leg.Draw();

        c.RedrawAxis();

        const auto out = plot_output_file(output_stem).string();
        c.SaveAs(out.c_str());
        std::cout << "\n[plot_roc_by_channel] saved plot: " << out << "\n";
        std::cout << "[plot_roc_by_channel] done\n";

        return 0;
    });
}

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif
