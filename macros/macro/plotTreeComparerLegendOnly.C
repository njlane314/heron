/* -- C++ -- */
/**
 * @file plot/macro/plotTreeComparerLegendOnly.C
 *
 * @brief Draw only the analysis-channel legend used by stacked histograms.
 *
 * Usage:
 *   root -l -q 'plot/macro/plotTreeComparerLegendOnly.C("analysis_channel_legend_only.pdf")'
 */

#include <memory>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"

#include "../include/PlotChannels.hh"
#include "MacroGuard.hh"

using namespace nu;

void plotTreeComparerLegendOnly(const char *output_name = "analysis_channel_legend_only.pdf")
{
  heron::macro::run_with_guard("plotTreeComparerLegendOnly", [&]() {
    gStyle->SetOptStat(0);

    TCanvas canvas("c_analysis_channel_legend", "Analysis channel legend only", 1200, 180);
    canvas.SetFillColor(kWhite);
    canvas.SetFrameFillColor(kWhite);

    TLegend legend(0.01, 0.05, 0.99, 0.95);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.14);
    legend.SetNColumns(4);
    legend.SetEntrySeparation(0.01);
    legend.SetColumnSeparation(0.04);

    std::vector<std::unique_ptr<TH1D>> proxies;
    const std::vector<int> channel_keys = Channels::mc_keys();
    proxies.reserve(channel_keys.size());

    for (size_t i = 0; i < channel_keys.size(); ++i)
    {
        const int channel_key = channel_keys[i];
        std::unique_ptr<TH1D> proxy(new TH1D(Form("proxy_%d", channel_key), "", 1, 0.0, 1.0));
        proxy->SetDirectory(nullptr);
        proxy->SetFillColor(Channels::colour(channel_key));
        proxy->SetFillStyle(Channels::fill_style(channel_key));
        proxy->SetLineColor(kBlack);
        proxy->SetLineWidth(1);
        legend.AddEntry(proxy.get(), Channels::label(channel_key).c_str(), "f");
        proxies.push_back(std::move(proxy));
    }

    legend.Draw();
    canvas.SaveAs(output_name);

  });
}
