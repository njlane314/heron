/* -- C++ -- */
/**
 * @file plot/macro/plotLegendOnlyNoStats.C
 *
 * @brief Draw only the stacked-category legend (no statistics box).
 *
 * Usage:
 *   root -l -q 'plot/macro/plotLegendOnlyNoStats.C("legend_only_no_stats.pdf")'
 */

#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "include/MacroGuard.hh"

namespace
{
struct LegendEntry
{
    std::string label;
    Color_t colour;
    int fill_style;
};

std::vector<LegendEntry> default_entries()
{
    std::vector<LegendEntry> entries;
    entries.push_back({"Out FV", TColor::GetColor("#ff9f1c"), 1001});
    entries.push_back({"#nu_{#mu}CC #pi^{0}/#gamma#gamma", TColor::GetColor("#ff006e"), 1001});
    entries.push_back({"#nu_{#mu}CC 0p1#pi^{#pm}", TColor::GetColor("#8338ec"), 1001});
    entries.push_back({"#nu_{#mu}CC Np0#pi", TColor::GetColor("#3a86ff"), 1001});
    entries.push_back({"#nu_{#mu}CC multi-#pi^{#pm}", TColor::GetColor("#5e60ce"), 1001});
    entries.push_back({"#nu_{x}NC", TColor::GetColor("#4361ee"), 1001});
    entries.push_back({"Cosmic", TColor::GetColor("#1f2a44"), 3345});
    entries.push_back({"#nu_{#mu}CC Other", TColor::GetColor("#118ab2"), 1001});
    LegendEntry dirt_entry = {"Dirt", kYellow, 1001};
    dirt_entry.colour = TColor::GetColor("#ffb703");
    entries.push_back(dirt_entry);
    entries.push_back({"Other", TColor::GetColor("#ff3d00"), 1001});
    entries.push_back({"Signal #Lambda^{0} CCQE (#Lambda^{0} #rightarrow p#pi^{-})", TColor::GetColor("#39ff14"), 1001});
    entries.push_back({"Signal #Lambda^{0} CCRES (#Lambda^{0} #rightarrow p#pi^{-})", TColor::GetColor("#00e676"), 1001});
    entries.push_back({"Signal #Lambda^{0} CCDIS (#Lambda^{0} #rightarrow p#pi^{-})", TColor::GetColor("#7ae582"), 1001});
    entries.push_back({"Signal #Lambda^{0} CC Other (#Lambda^{0} #rightarrow p#pi^{-})", TColor::GetColor("#2dc653"), 1001});
    return entries;
}
} // namespace

void plotLegendOnlyNoStats(const char *output_name = "legend_only_no_stats.pdf")
{
  heron::macro::run_with_guard("plotLegendOnlyNoStats", [&]() {
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

  });
}
