/* -- C++ -- */
/**
 * @file macros/plot_legend.C
 *
 * @brief Draw only the stacked-category legend (no statistics box).
 *
 * Usage:
 *   root -l -q 'plot/macro/plot_legend.C("plot_legend.pdf")'
 */

#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"

#include "AnalysisChannels.hh"
#include "PlotChannels.hh"

namespace
{
struct LegendEntry
{
    std::string label;
    Color_t colour;
    int fill_style;
};

std::vector<int> channel_codes()
{
    return {
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::OutFV),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::MuCCPi0OrGamma),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::MuCC1pi),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::MuCC0pi_ge1p),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::MuCCNpi),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::NC),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::External),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::MuCCOther),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::MuCCSigma0),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::MuCCK0),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::ECCC),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::SignalLambda),
        AnalysisChannels::to_int(AnalysisChannels::AnalysisChannel::DataInclusive)
    };
}

std::vector<LegendEntry> default_entries()
{
    std::vector<LegendEntry> entries;
    const std::vector<int> channels = channel_codes();
    entries.reserve(channels.size());

    for (const int channel : channels)
    {
        const nu::Channels::Properties &properties = nu::Channels::properties(channel);
        entries.push_back({properties.tex_label, properties.fill_colour, properties.fill_style});
    }

    return entries;
}
} // namespace

void plot_legend(const char *output_name = "plot_legend.pdf")
{
    gStyle->SetOptStat(0);

    TCanvas canvas("c_legend_only", "Legend only", 700, 360);
    canvas.SetFillColor(kWhite);
    canvas.SetFrameFillColor(kWhite);

    TLegend legend(0.01, 0.05, 0.99, 0.95);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetLineColorAlpha(kWhite, 0.0);
    legend.SetLineWidth(0);
    legend.SetShadowColor(kWhite);
    legend.SetTextFont(42);
    legend.SetTextSize(0.08);
    legend.SetNColumns(2);
    legend.SetEntrySeparation(0.01);
    legend.SetColumnSeparation(0.04);

    const std::vector<LegendEntry> entries = default_entries();
    std::vector<TH1D *> proxies;
    proxies.reserve(entries.size());

    for (size_t i = 0; i < entries.size(); ++i)
    {
        TH1D *proxy = new TH1D(Form("proxy_%zu", i), "", 1, 0.0, 1.0);
        proxy->SetDirectory(nullptr);
        proxy->SetFillColor(entries[i].colour);
        proxy->SetFillStyle(entries[i].fill_style);
        proxy->SetLineColor(kBlack);
        proxy->SetLineWidth(1);
        proxies.push_back(proxy);
        legend.AddEntry(proxy, entries[i].label.c_str(), "f");
    }

    legend.Draw();
    canvas.SaveAs(output_name);

    for (TH1D *proxy : proxies)
    {
        delete proxy;
    }
}
