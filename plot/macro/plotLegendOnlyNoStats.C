/* -- C++ -- */
/**
 * @file plot/macro/plotLegendOnlyNoStats.C
 *
 * @brief Draw only the stacked-category legend (no statistics box).
 *
 * Usage:
 *   root -l -q 'plot/macro/plotLegendOnlyNoStats.C("legend_only.png")'
 */

#include <string>
#include <vector>

#include "TH1D.h"
#include "TCanvas.h"
#include "TColor.h"
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
    entries.push_back({"Out FV : 20,859.04", kYellow - 7, 1001});
    entries.push_back({"#nu_{#mu}CC #pi^{0}/#gamma#gamma : 1,781.26", kOrange, 1001});
    entries.push_back({"#nu_{#mu}CC 0p1#pi^{#pm} : 5,680.95", kRed - 7, 1001});
    entries.push_back({"#nu_{#mu}CC Np0#pi : 4,486.91", kRed, 1001});
    entries.push_back({"#nu_{#mu}CC multi-#pi^{#pm} : 328.67", kViolet, 1001});
    entries.push_back({"#nu_{x}NC : 3,536.21", kBlue, 1001});
    entries.push_back({"#nu_{#mu}CC multi-strange : 73.41", kGreen + 2, 1001});
    entries.push_back({"Cosmic : 0.00", kTeal + 2, 3345});
    entries.push_back({"#nu_{#mu}CC Other : 468.59", kCyan + 2, 1001});
    entries.push_back({"Dirt : 1,111.52", TColor::GetColor("#f6d32d"), 1001});
    entries.push_back({"Other : 0.00", kCyan, 1001});
    entries.push_back({"#nu_{#mu}CC single-strange : 83.05", kSpring + 5, 1001});
    return entries;
}
} // namespace

void plotLegendOnlyNoStats(const char *output_name = "legend_only_no_stats.png")
{
    gStyle->SetOptStat(0);

    TCanvas canvas("c_legend_only", "Legend only", 1250, 220);
    canvas.SetFillColor(kWhite);
    canvas.SetFrameFillColor(kWhite);

    TH1D frame("frame", "", 1, 0.0, 1.0);
    frame.SetStats(0);
    frame.SetDirectory(nullptr);
    frame.GetXaxis()->SetLabelSize(0.0);
    frame.GetYaxis()->SetLabelSize(0.0);
    frame.GetXaxis()->SetTickLength(0.0);
    frame.GetYaxis()->SetTickLength(0.0);
    frame.SetLineColor(kWhite);
    frame.Draw();

    TLegend legend(0.01, 0.05, 0.99, 0.95);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
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
