/* -- C++ -- */
/// \file plot/macro/plotPotSimple.C
/// \brief POT timeline plotting macro for Nuxsec.

#include "RVersion.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include <sqlite3.h>

#include "../include/Plotter.hh"
#include "../include/PlotEnv.hh"

namespace
{

std::string db_root()
{
    if (const char *e = gSystem->Getenv("DBROOT"))
    {
        return e;
    }
    return "/exp/uboone/data/uboonebeam/beamdb";
}

std::string slip_dir()
{
    if (const char *e = gSystem->Getenv("SLIP_DIR"))
    {
        return e;
    }
    return "/exp/uboone/app/users/guzowski/slip_stacking";
}

time_t iso_to_utc(const char *s)
{
    std::tm tm{};
    strptime(s, "%Y-%m-%dT%H:%M:%S", &tm);
    return timegm(&tm);
}

time_t sunday_after_or_on(time_t t)
{
    std::tm gm = *gmtime(&t);
    const int add = (7 - gm.tm_wday) % 7;
    gm.tm_hour = gm.tm_min = gm.tm_sec = 0;
    const time_t day0 = timegm(&gm);
    return day0 + add * 86400;
}

void configure_style()
{
    nuxsec::plot::Plotter{}.set_global_style();
    // We draw a custom (blue) right axis; avoid black RHS ticks.
    gStyle->SetPadTickY(0);
    TGaxis::SetMaxDigits(3);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
    gStyle->SetTimeOffset(0, "gmt");
#else
    gStyle->SetTimeOffset(0);
#endif
}

sqlite3 *open_database(const std::string &run_db)
{
    sqlite3 *db = nullptr;
    if (sqlite3_open(run_db.c_str(), &db) != SQLITE_OK)
    {
        std::cerr << "Could not open " << run_db << "\n";
        return nullptr;
    }
    return db;
}

void attach_database(sqlite3 *db, const std::string &path, const std::string &alias)
{
    const auto sql = "ATTACH DATABASE '" + path + "' AS " + alias + ";";
    char *err = nullptr;
    const int rc = sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &err);
    if (rc != SQLITE_OK && err)
    {
        std::cerr << "SQL error: " << err << "\n";
        sqlite3_free(err);
    }
}

struct pot_samples
{
    std::vector<double> times;
    std::vector<double> pots;
};

pot_samples fetch_samples(sqlite3 *db, const char *sql)
{
    pot_samples samples;
    sqlite3_stmt *statement = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &statement, nullptr) != SQLITE_OK)
    {
        return samples;
    }
    while (sqlite3_step(statement) == SQLITE_ROW)
    {
        const char *begin_time = reinterpret_cast<const char *>(sqlite3_column_text(statement, 0));
        const double pot = sqlite3_column_double(statement, 1);
        if (!begin_time || pot <= 0)
        {
            continue;
        }
        samples.times.push_back(static_cast<double>(iso_to_utc(begin_time)));
        samples.pots.push_back(pot);
    }
    sqlite3_finalize(statement);
    return samples;
}

bool determine_range(const std::vector<pot_samples> &samples, double &xlo, double &xhi, int &nbins)
{
    double tmin = 0;
    double tmax = 0;
    bool have_range = false;
    for (const auto &sample : samples)
    {
        for (double t : sample.times)
        {
            if (!have_range)
            {
                tmin = t;
                tmax = t;
                have_range = true;
            }
            else
            {
                tmin = std::min(tmin, t);
                tmax = std::max(tmax, t);
            }
        }
    }
    if (!have_range)
    {
        return false;
    }
    const double week = 7.0 * 86400.0;
    const time_t first_sunday = sunday_after_or_on(static_cast<time_t>(tmin));
    const time_t last_sunday = sunday_after_or_on(static_cast<time_t>(tmax));
    xlo = static_cast<double>(first_sunday - 7 * 86400);
    xhi = static_cast<double>(last_sunday + 7 * 86400);
    nbins = std::max(1, static_cast<int>((xhi - xlo) / week + 0.5));
    return true;
}

struct histogram_bundle
{
    TH1D bnb;
    TH1D fhc;
    TH1D rhc;

    histogram_bundle(int nbins, double xlo, double xhi)
        : bnb("hBNB", "", nbins, xlo, xhi)
        , fhc("hFHC", "", nbins, xlo, xhi)
        , rhc("hRHC", "", nbins, xlo, xhi)
    {
        bnb.SetDirectory(nullptr);
        fhc.SetDirectory(nullptr);
        rhc.SetDirectory(nullptr);

        // Brighter, bolder fills for the timeline.
        const Int_t col_bnb = TColor::GetColor("#00c853");
        const Int_t col_fhc = TColor::GetColor("#ffb300");
        const Int_t col_rhc = TColor::GetColor("#ff1744");

        bnb.SetFillColor(col_bnb);
        bnb.SetLineColor(kBlack);
        bnb.SetLineWidth(1);
        bnb.SetBarWidth(1.0);
        bnb.SetBarOffset(0.0);

        fhc.SetFillColor(col_fhc);
        fhc.SetLineColor(kBlack);
        fhc.SetLineWidth(1);
        fhc.SetBarWidth(1.0);
        fhc.SetBarOffset(0.0);

        rhc.SetFillColor(col_rhc);
        rhc.SetLineColor(kBlack);
        rhc.SetLineWidth(1);
        rhc.SetBarWidth(1.0);
        rhc.SetBarOffset(0.0);
    }
};

void fill_histogram(const pot_samples &samples, TH1D &histogram)
{
    for (size_t i = 0; i < samples.times.size(); ++i)
    {
        histogram.Fill(samples.times[i], samples.pots[i] / 1e18);
    }
}

struct cumulative_data
{
    std::vector<double> x;
    std::vector<double> scaled_total;
    std::vector<double> scaled_bnb;
    std::vector<double> scaled_fhc;
    std::vector<double> scaled_rhc;
    double y_max;
    double max_cumulative_total;
};

cumulative_data compute_cumulative_data(const histogram_bundle &histograms, int nbins)
{
    cumulative_data data;
    data.x.resize(nbins);
    data.scaled_total.resize(nbins);
    data.scaled_bnb.resize(nbins);
    data.scaled_fhc.resize(nbins);
    data.scaled_rhc.resize(nbins);

    std::vector<double> cum_total(nbins);
    std::vector<double> cum_bnb(nbins);
    std::vector<double> cum_fhc(nbins);
    std::vector<double> cum_rhc(nbins);

    double max_stack = 0;
    double sum_total_pot = 0;
    double sum_bnb_pot = 0;
    double sum_fhc_pot = 0;
    double sum_rhc_pot = 0;
    data.max_cumulative_total = 0;
    for (int i = 1; i <= nbins; ++i)
    {
        const double w_bnb = histograms.bnb.GetBinContent(i);
        const double w_fhc = histograms.fhc.GetBinContent(i);
        const double w_rhc = histograms.rhc.GetBinContent(i);
        const double stack = w_bnb + w_fhc + w_rhc;
        max_stack = std::max(max_stack, stack);

        // Histogram bin contents are in units of 1e18 POT/week; integrate in POT then scale to 1e20.
        sum_total_pot += stack * 1e18;
        sum_bnb_pot += w_bnb * 1e18;
        sum_fhc_pot += w_fhc * 1e18;
        sum_rhc_pot += w_rhc * 1e18;

        const double c_total = sum_total_pot / 1e20;
        const double c_bnb = sum_bnb_pot / 1e20;
        const double c_fhc = sum_fhc_pot / 1e20;
        const double c_rhc = sum_rhc_pot / 1e20;

        cum_total[i - 1] = c_total;
        cum_bnb[i - 1] = c_bnb;
        cum_fhc[i - 1] = c_fhc;
        cum_rhc[i - 1] = c_rhc;

        data.max_cumulative_total = std::max(data.max_cumulative_total, c_total);
        data.x[i - 1] = histograms.bnb.GetXaxis()->GetBinCenter(i);
    }
    // "Nice" headroom like the reference (major ticks at 0,5,10,... with minor ticks between).
    if (max_stack > 0)
    {
        const double nice = std::ceil(max_stack / 5.0) * 5.0;
        data.y_max = nice * 1.08;
    }
    else data.y_max = 1.0;
    // Scale all cumulative curves to the left-axis range using the TOTAL cumulative maximum.
    const double scale = data.max_cumulative_total > 0 ? data.y_max / data.max_cumulative_total : 1.0;
    for (int i = 0; i < nbins; ++i)
    {
        data.scaled_total[i] = cum_total[i] * scale;
        data.scaled_bnb[i] = cum_bnb[i] * scale;
        data.scaled_fhc[i] = cum_fhc[i] * scale;
        data.scaled_rhc[i] = cum_rhc[i] * scale;
    }
    return data;
}

void draw_plot(const histogram_bundle &histograms, const cumulative_data &data, const std::string &outfile)
{
    // Aspect and margins closer to the reference figure.
    TCanvas canvas("c", "POT timeline",
                   nuxsec::plot::kCanvasWidth,
                   nuxsec::plot::kCanvasHeight);
    const double ml = 0.08;
    const double mr = 0.12;
    const double split = 0.82;
    TPad *main_pad = new TPad("pad_main", "pad_main", 0.0, 0.0, 1.0, split);
    TPad *legend_pad = new TPad("pad_legend", "pad_legend", 0.0, split, 1.0, 1.0);
    main_pad->SetTopMargin(0.02);
    main_pad->SetBottomMargin(0.20);
    main_pad->SetLeftMargin(ml);
    main_pad->SetRightMargin(mr);
    legend_pad->SetTopMargin(0.08);
    legend_pad->SetBottomMargin(0.02);
    // Match main-pad margins so the legend aligns with the plot frame (between y-axes).
    legend_pad->SetLeftMargin(ml);
    legend_pad->SetRightMargin(mr);
    main_pad->Draw();
    legend_pad->Draw();
    main_pad->cd();
    main_pad->SetGridy(false);
    main_pad->SetTickx(1);
    main_pad->SetTicky(0);

    THStack stack("hs", "");
    stack.Add(const_cast<TH1D *>(&histograms.bnb));
    stack.Add(const_cast<TH1D *>(&histograms.fhc));
    stack.Add(const_cast<TH1D *>(&histograms.rhc));
    stack.Draw("hist");
    stack.GetXaxis()->SetTimeDisplay(1);
    stack.GetXaxis()->SetTimeOffset(0, "gmt");
    stack.GetXaxis()->SetTimeFormat("%d/%b/%Y");
    stack.GetXaxis()->SetNdivisions(509);
    stack.GetXaxis()->SetLabelSize(0.045);
    stack.GetXaxis()->SetLabelOffset(0.018);

    stack.GetYaxis()->SetTitle("POT per week (x 10^{18})");
    stack.GetYaxis()->SetNdivisions(507);
    stack.GetYaxis()->SetTitleSize(0.040);
    stack.GetYaxis()->SetLabelSize(0.045);
    stack.GetYaxis()->SetTitleOffset(0.70);

    stack.SetMaximum(data.y_max);
    stack.SetMinimum(0);

    // Cumulative curves (scaled onto the left axis).
    const Int_t col_total = TColor::GetColor("#2196f3");
    const Int_t col_bnb = TColor::GetColor("#00e676");
    const Int_t col_fhc = TColor::GetColor("#ffc400");
    const Int_t col_rhc = TColor::GetColor("#ff1744");

    auto style_line = [](TGraph &g, Int_t col, Int_t lstyle, Int_t lwidth,
                         Int_t mstyle, double msize, bool markers)
    {
        g.SetLineColor(col);
        g.SetLineStyle(lstyle);
        g.SetLineWidth(lwidth);
        if (markers)
        {
            g.SetMarkerColor(col);
            g.SetMarkerStyle(mstyle);
            g.SetMarkerSize(msize);
        }
    };

    // Total cumulative POT: keep it dominant and clean (solid blue) with an outline for contrast.
    TGraph g_total(data.x.size(), data.x.data(), data.scaled_total.data());
    TGraph g_total_outline(data.x.size(), data.x.data(), data.scaled_total.data());
    style_line(g_total_outline, kBlack, 1, 5, 20, 0.0, false);
    style_line(g_total, col_total, 1, 4, 20, 0.0, false);
    g_total_outline.Draw("L SAME");
    g_total.Draw("L SAME");

    // Beam-mode cumulative curves: thicker, distinct dash patterns + markers + outline.
    // (Markers help against filled histograms; outlines avoid low-contrast segments.)
    TGraph g_bnb(data.x.size(), data.x.data(), data.scaled_bnb.data());
    TGraph g_bnb_outline(data.x.size(), data.x.data(), data.scaled_bnb.data());
    style_line(g_bnb_outline, kBlack, 1, 5, 21, 0.0, false);
    style_line(g_bnb, col_bnb, 1, 4, 21, 0.55, true);
    g_bnb_outline.Draw("L SAME");
    g_bnb.Draw("LP SAME");

    TGraph g_fhc(data.x.size(), data.x.data(), data.scaled_fhc.data());
    TGraph g_fhc_outline(data.x.size(), data.x.data(), data.scaled_fhc.data());
    style_line(g_fhc_outline, kBlack, 1, 5, 22, 0.0, false);
    style_line(g_fhc, col_fhc, 1, 4, 22, 0.55, true);
    g_fhc_outline.Draw("L SAME");
    g_fhc.Draw("LP SAME");

    TGraph g_rhc(data.x.size(), data.x.data(), data.scaled_rhc.data());
    TGraph g_rhc_outline(data.x.size(), data.x.data(), data.scaled_rhc.data());
    style_line(g_rhc_outline, kBlack, 1, 5, 23, 0.0, false);
    style_line(g_rhc, col_rhc, 1, 4, 23, 0.55, true);
    g_rhc_outline.Draw("L SAME");
    g_rhc.Draw("LP SAME");

    const double xhi = stack.GetXaxis()->GetXmax();
    const double rhs_max = (data.max_cumulative_total > 0) ? (data.max_cumulative_total * 1.03) : 1.0;
    TGaxis right_axis(xhi, 0, xhi, data.y_max, 0, rhs_max, 507, "+L");
    right_axis.SetLineColor(col_total);
    right_axis.SetLabelColor(col_total);
    right_axis.SetTitleColor(col_total);
    right_axis.SetLineWidth(2);
    right_axis.SetLabelFont(42);
    right_axis.SetTitleFont(42);
    right_axis.SetLabelSize(0.045);
    right_axis.SetTitleSize(0.040);
    right_axis.SetTitleOffset(0.80);
    right_axis.SetTickSize(0.02);
    right_axis.SetTitle("Cumulative POT (x 10^{20})");
    right_axis.Draw();

    legend_pad->cd();
    // Keep legend strictly within the plot-frame x-extent (between LHS and RHS axes).
    const double lx1 = legend_pad->GetLeftMargin() + 0.02;
    const double lx2 = 1.0 - legend_pad->GetRightMargin() - 0.02;
    const double ly1 = legend_pad->GetBottomMargin() + 0.02;
    const double ly2 = 1.0 - legend_pad->GetTopMargin() - 0.02;
    TLegend legend(lx1, ly1, lx2, ly2);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetNColumns(2);
    legend.SetColumnSeparation(0.08);
    legend.SetEntrySeparation(0.00);
    legend.SetMargin(0.25);
    const double s_main = 0.045;
    const double s_leg = s_main * (split / (1.0 - split));
    legend.SetTextSize(s_leg);
    legend.AddEntry(&histograms.bnb, "BNB (\\nu)", "f");
    legend.AddEntry(&histograms.fhc, "NuMI-FHC (\\nu)", "f");
    legend.AddEntry(&histograms.rhc, "NuMI-RHC (\\bar{\\nu})", "f");
    legend.AddEntry(&g_total, "Total cumulative POT", "l");
    legend.AddEntry(&g_bnb, "BNB cumulative POT", "lp");
    legend.AddEntry(&g_fhc, "NuMI-FHC cumulative POT", "lp");
    legend.AddEntry(&g_rhc, "NuMI-RHC cumulative POT", "lp");
    legend.Draw();
    main_pad->cd();
    main_pad->RedrawAxis();
    canvas.SaveAs(outfile.c_str());
    delete main_pad;
    delete legend_pad;
}

void plotPotSimpleInternal(const char *out = nullptr)
{
    configure_style();
    const std::string run_db = db_root() + "/run.db";
    const std::string bnb_db = gSystem->AccessPathName((db_root() + "/bnb_v2.db").c_str())
                                   ? db_root() + "/bnb_v1.db"
                                   : db_root() + "/bnb_v2.db";
    const std::string numi_db = gSystem->AccessPathName((db_root() + "/numi_v2.db").c_str())
                                    ? db_root() + "/numi_v1.db"
                                    : db_root() + "/numi_v2.db";
    const std::string n4_db = slip_dir() + "/numi_v4.db";
    sqlite3 *db = open_database(run_db);
    if (!db)
    {
        return;
    }
    attach_database(db, bnb_db, "bnb");
    attach_database(db, numi_db, "numi");
    attach_database(db, n4_db, "n4");
    const char *q_bnb =
        "SELECT r.begin_time, 1e12 * ("
        " CASE WHEN IFNULL(b.tor875,0)>0 THEN IFNULL(b.tor875,r.tor875) "
        "      WHEN IFNULL(r.tor875,0)>0 THEN r.tor875 "
        "      WHEN IFNULL(b.tor860,0)>0 THEN IFNULL(b.tor860,r.tor860) "
        "      ELSE IFNULL(r.tor860,0) END ) AS pot "
        "FROM runinfo r LEFT JOIN bnb.bnb b ON r.run=b.run AND r.subrun=b.subrun "
        "WHERE (IFNULL(b.tor875,0)+IFNULL(r.tor875,0)+IFNULL(b.tor860,0)+IFNULL(r.tor860,0))>0;";
    const char *q_fhc =
        "SELECT r.begin_time, 1e12 * (CASE WHEN IFNULL(n.tortgt_fhc,0)>0 THEN n.tortgt_fhc "
        "                                  ELSE IFNULL(n.tor101_fhc,0) END) AS pot "
        "FROM runinfo r JOIN n4.numi n ON r.run=n.run AND r.subrun=n.subrun "
        "WHERE IFNULL(n.EA9CNT_fhc,0)>0;";
    const char *q_rhc =
        "SELECT r.begin_time, 1e12 * (CASE WHEN IFNULL(n.tortgt_rhc,0)>0 THEN n.tortgt_rhc "
        "                                  ELSE IFNULL(n.tor101_rhc,0) END) AS pot "
        "FROM runinfo r JOIN n4.numi n ON r.run=n.run AND r.subrun=n.subrun "
        "WHERE IFNULL(n.EA9CNT_rhc,0)>0;";
    const pot_samples bnb_samples = fetch_samples(db, q_bnb);
    const pot_samples fhc_samples = fetch_samples(db, q_fhc);
    const pot_samples rhc_samples = fetch_samples(db, q_rhc);
    sqlite3_close(db);
    double xlo = 0;
    double xhi = 0;
    int nbins = 0;
    const std::vector<pot_samples> all_samples{bnb_samples, fhc_samples, rhc_samples};
    if (!determine_range(all_samples, xlo, xhi, nbins))
    {
        std::cerr << "No time range\n";
        return;
    }
    histogram_bundle histograms(nbins, xlo, xhi);
    fill_histogram(bnb_samples, histograms.bnb);
    fill_histogram(fhc_samples, histograms.fhc);
    fill_histogram(rhc_samples, histograms.rhc);
    const cumulative_data data = compute_cumulative_data(histograms, nbins);
    const auto out_path = nuxsec::plot::resolve_output_file(out ? out : "", "pot_timeline");
    const std::string outfile = out_path.string();
    draw_plot(histograms, data, outfile);
}

} // namespace

int nuxsec_plot(const char *out = nullptr)
{
    plotPotSimpleInternal(out);
    return 0;
}

void plotPotSimple(const char *out = nullptr)
{
    plotPotSimpleInternal(out);
}
