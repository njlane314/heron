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
    entries.push_back({"Out FV", kYellow - 7, 1001});
    entries.push_back({"#nu_{#mu}CC #pi^{0}/#gamma#gamma", kOrange, 1001});
    entries.push_back({"#nu_{#mu}CC 0p1#pi^{#pm}", kRed - 7, 1001});
    entries.push_back({"#nu_{#mu}CC Np0#pi", kRed, 1001});
    entries.push_back({"#nu_{#mu}CC multi-#pi^{#pm}", kViolet, 1001});
    entries.push_back({"#nu_{x}NC", kBlue, 1001});
    entries.push_back({"#nu_{#mu}CC multi-strange", kGreen + 2, 1001});
    entries.push_back({"Cosmic", kTeal + 2, 3345});
    entries.push_back({"#nu_{#mu}CC Other", kCyan + 2, 1001});
    LegendEntry dirt_entry = {"Dirt", kYellow, 1001};
    dirt_entry.colour = TColor::GetColor("#f6d32d");
    entries.push_back(dirt_entry);
    entries.push_back({"Other", kCyan, 1001});
    entries.push_back({"#nu_{#mu}CC single-strange", kSpring + 5, 1001});
    return entries;
}
} // namespace

void plotLegendOnlyNoStats(const char *output_name = "legend_only_no_stats.pdf")
{
    gStyle->SetOptStat(0);

    TCanvas canvas("c_legend_only", "Legend only", 600, 120);
    canvas.SetFillColor(kWhite);
    canvas.SetFrameFillColor(kWhite);

    TLegend legend(0.01, 0.05, 0.99, 0.95);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetLineColorAlpha(kWhite, 0.0);
    legend.SetLineWidth(0);
    legend.SetShadowColor(kWhite);
    legend.SetTextFont(42);
    legend.SetTextSize(0.14);
    legend.SetNColumns(4);
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
