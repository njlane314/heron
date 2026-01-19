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
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(1);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelSize(0.038, "XYZ");
    gStyle->SetTitleSize(0.050, "YZ");
    gStyle->SetTitleOffset(1.00, "Y");
    gStyle->SetPadTickX(1);
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
    std::vector<double> scaled;
    double y_max;
    double max_cumulative;
};

cumulative_data compute_cumulative_data(const histogram_bundle &histograms, int nbins)
{
    cumulative_data data;
    data.x.resize(nbins);
    data.scaled.resize(nbins);
    std::vector<double> cumulative(nbins);
    double max_stack = 0;
    double sum = 0;
    data.max_cumulative = 0;
    for (int i = 1; i <= nbins; ++i)
    {
        const double stack = histograms.bnb.GetBinContent(i) + histograms.fhc.GetBinContent(i) +
                             histograms.rhc.GetBinContent(i);
        max_stack = std::max(max_stack, stack);
        sum += stack * 1e18;
        const double cumulative_value = sum / 1e20;
        cumulative[i - 1] = cumulative_value;
        data.max_cumulative = std::max(data.max_cumulative, cumulative_value);
        data.x[i - 1] = histograms.bnb.GetXaxis()->GetBinCenter(i);
    }
    // "Nice" headroom like the reference (major ticks at 0,5,10,... with minor ticks between).
    if (max_stack > 0)
    {
        const double nice = std::ceil(max_stack / 5.0) * 5.0;
        data.y_max = nice * 1.02;
    }
    else data.y_max = 1.0;
    const double scale = data.max_cumulative > 0 ? data.y_max / data.max_cumulative : 1.0;
    for (int i = 0; i < nbins; ++i)
    {
        data.scaled[i] = cumulative[i] * scale;
    }
    return data;
}

void draw_plot(const histogram_bundle &histograms, const cumulative_data &data, const std::string &outfile)
{
    // Aspect and margins closer to the reference figure.
    TCanvas canvas("c", "POT timeline", 1200, 700);
    canvas.SetMargin(0.08, 0.12, 0.20, 0.08);
    canvas.SetGridy(false);
    canvas.SetTickx(1);
    canvas.SetTicky(0);

    THStack stack("hs", "");
    stack.Add(const_cast<TH1D *>(&histograms.bnb));
    stack.Add(const_cast<TH1D *>(&histograms.fhc));
    stack.Add(const_cast<TH1D *>(&histograms.rhc));
    stack.Draw("hist");
    stack.GetXaxis()->SetTimeDisplay(1);
    stack.GetXaxis()->SetTimeOffset(0, "gmt");
    stack.GetXaxis()->SetTimeFormat("%d/%b/%Y");
    stack.GetXaxis()->SetNdivisions(509);
    stack.GetXaxis()->SetLabelSize(0.032);
    stack.GetXaxis()->SetLabelOffset(0.018);

    stack.GetYaxis()->SetTitle("POT per week (x 10^{18})");
    stack.GetYaxis()->SetNdivisions(507);
    stack.GetYaxis()->SetTitleSize(0.050);
    stack.GetYaxis()->SetLabelSize(0.036);
    stack.GetYaxis()->SetTitleOffset(0.70);

    stack.SetMaximum(data.y_max);
    stack.SetMinimum(0);

    const Int_t col_cumulative = TColor::GetColor("#2962ff");
    TGraph graph(data.x.size(), data.x.data(), data.scaled.data());
    graph.SetLineColor(col_cumulative);
    graph.SetLineWidth(2);
    graph.Draw("L SAME");

    const double xhi = stack.GetXaxis()->GetXmax();
    TGaxis right_axis(xhi, 0, xhi, data.y_max, 0, data.max_cumulative, 507, "+L");
    right_axis.SetLineColor(col_cumulative);
    right_axis.SetLabelColor(col_cumulative);
    right_axis.SetTitleColor(col_cumulative);
    right_axis.SetLineWidth(2);
    right_axis.SetLabelFont(42);
    right_axis.SetTitleFont(42);
    right_axis.SetLabelSize(0.036);
    right_axis.SetTitleSize(0.050);
    right_axis.SetTitleOffset(0.80);
    right_axis.SetTickSize(0.02);
    right_axis.SetTitle("Cumulative POT (x 10^{20})");
    right_axis.Draw();

    TLegend legend(0.16, 0.68, 0.40, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.034);
    legend.AddEntry(&histograms.bnb, "BNB (\\nu)", "f");
    legend.AddEntry(&histograms.fhc, "NuMI-FHC (\\nu)", "f");
    legend.AddEntry(&histograms.rhc, "NuMI-RHC (\\bar{\\nu})", "f");
    legend.AddEntry(&graph, "Total POT", "l");
    legend.Draw();
    canvas.RedrawAxis();
    canvas.SaveAs(outfile.c_str());
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
