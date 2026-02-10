/* -- C++ -- */
/// \file plot/macro/plotPotTimelineStyle.C
/// \brief MicroBooNE POT timeline with cumulative POT, beam power, and run bands:
///        cumulative POT (left axis) + beam power (right axis) + run bands.
///
/// Notes:
///  - Uses the same MicroBooNE sqlite beam DBs + queries as plotPotSimple.C.
///  - Beam power points are computed as *average* power over the time bin:
///      P[kW] = POT/bin_width[s] * E[GeV] * 1.602176634e-10[J/GeV] / 1000
///  - Binning defaults to 1 week (Sunday-to-Sunday, UTC) to match your existing macro.
///
/// Usage:
///   root -l -q 'plotPotTimelineStyle.C("pot_timeline_style.pdf")'
/// or within your framework:
///   root -l -q 'plotPotTimelineStyle.C'
///
/// Env:
///   DBROOT   (optional) e.g. /exp/uboone/data/uboonebeam/beamdb
///   SLIP_DIR (optional) e.g. /exp/uboone/app/users/guzowski/slip_stacking

#include "RVersion.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <sqlite3.h>

#include "../include/Plotter.hh"
#include "../include/PlotEnv.hh"

using namespace nu;


// -------------------- Paths --------------------

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

// -------------------- Time helpers --------------------

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

std::string utc_date_ymd(time_t t)
{
    char buf[64];
    std::tm gm = *gmtime(&t);
    strftime(buf, sizeof(buf), "%Y-%m-%d", &gm);
    return std::string(buf);
}

// -------------------- Style --------------------

void configure_style()
{
    Plotter{}.set_global_style();
    // We'll draw a custom right axis; avoid default RHS ticks.
    gStyle->SetPadTickY(0);
    // Subtle dotted grid (closer to typical "timeline" figures).
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(kGray + 1);
    gStyle->SetGridWidth(1);
    TGaxis::SetMaxDigits(3);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
    gStyle->SetTimeOffset(0, "gmt");
#else
    gStyle->SetTimeOffset(0);
#endif
}

// -------------------- SQLite helpers --------------------

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

// -------------------- Data fetch --------------------

struct pot_samples
{
    std::vector<double> times; // unix seconds (UTC)
    std::vector<double> pots;  // POT (protons)
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

bool determine_range_weekly(const std::vector<pot_samples> &samples, double &xlo, double &xhi, int &nbins,
                            time_t &tmin_out, time_t &tmax_out)
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

    // Weekly bins (Sunday-to-Sunday, UTC), padded by a week on both sides.
    const double week = 7.0 * 86400.0;
    const time_t first_sunday = sunday_after_or_on(static_cast<time_t>(tmin));
    const time_t last_sunday = sunday_after_or_on(static_cast<time_t>(tmax));

    xlo = static_cast<double>(first_sunday - 7 * 86400);
    xhi = static_cast<double>(last_sunday + 7 * 86400);

    nbins = std::max(1, static_cast<int>((xhi - xlo) / week + 0.5));

    tmin_out = static_cast<time_t>(tmin);
    tmax_out = static_cast<time_t>(tmax);
    return true;
}

void fill_histogram_pot(const pot_samples &samples, TH1D &histogram)
{
    for (size_t i = 0; i < samples.times.size(); ++i)
    {
        histogram.Fill(samples.times[i], samples.pots[i]);
    }
}

// -------------------- Run-period inference (from binned total POT) --------------------

struct run_period
{
    double x1;          // start (unix seconds)
    double x2;          // end   (unix seconds)
    std::string label;  // e.g. "Run1"
};

std::vector<run_period> infer_run_periods(const TH1D &h_bnb, const TH1D &h_fhc, const TH1D &h_rhc,
                                         double active_threshold_pot,
                                         double merge_gap_seconds)
{
    std::vector<std::pair<int, int>> segments; // [start_bin, end_bin], inclusive
    const int nbins = h_bnb.GetNbinsX();

    bool in_seg = false;
    int start = -1;

    for (int i = 1; i <= nbins; ++i)
    {
        const double pot = h_bnb.GetBinContent(i) + h_fhc.GetBinContent(i) + h_rhc.GetBinContent(i);
        const bool active = (pot > active_threshold_pot);

        if (active && !in_seg)
        {
            in_seg = true;
            start = i;
        }
        if (!active && in_seg)
        {
            in_seg = false;
            const int end = i - 1;
            if (end >= start)
            {
                segments.emplace_back(start, end);
            }
            start = -1;
        }
    }
    if (in_seg && start >= 1)
    {
        segments.emplace_back(start, nbins);
    }

    // Merge nearby segments separated by short inactive gaps.
    std::vector<std::pair<int, int>> merged;
    for (const auto &seg : segments)
    {
        if (merged.empty())
        {
            merged.push_back(seg);
            continue;
        }
        const auto &prev = merged.back();

        const double prev_x2 = h_bnb.GetXaxis()->GetBinUpEdge(prev.second);
        const double cur_x1 = h_bnb.GetXaxis()->GetBinLowEdge(seg.first);
        const double gap = cur_x1 - prev_x2;

        if (gap <= merge_gap_seconds)
        {
            merged.back().second = seg.second;
        }
        else
        {
            merged.push_back(seg);
        }
    }

    std::vector<run_period> out;
    out.reserve(merged.size());
    for (size_t i = 0; i < merged.size(); ++i)
    {
        const int b1 = merged[i].first;
        const int b2 = merged[i].second;
        run_period rp;
        rp.x1 = h_bnb.GetXaxis()->GetBinLowEdge(b1);
        rp.x2 = h_bnb.GetXaxis()->GetBinUpEdge(b2);
        rp.label = "Run" + std::to_string(i + 1);
        out.push_back(rp);
    }
    return out;
}

// -------------------- Plot data assembly --------------------

struct plot_data
{
    // Step curves: x has length (1 + 2*nbins)
    std::vector<double> x_step;
    std::vector<double> y_total;
    std::vector<double> y_bnb;
    std::vector<double> y_fhc;
    std::vector<double> y_rhc;

    // Beam power points (mapped onto left-axis coordinates; right axis is logarithmic)
    std::vector<double> x_pow_bnb;
    std::vector<double> y_pow_bnb_scaled;
    std::vector<double> x_pow_fhc;
    std::vector<double> y_pow_fhc_scaled;
    std::vector<double> x_pow_rhc;
    std::vector<double> y_pow_rhc_scaled;

    // Axis maxima
    double y_pot_max = 1.0;     // left axis (x 1e20)
    double pow_axis_min = 1.0;  // right axis (kW) - log scale (must be > 0)
    double pow_axis_max = 1.0;  // right axis (kW) - log scale

    // Totals (POT)
    double pot_total = 0;
    double pot_bnb = 0;
    double pot_fhc = 0;
    double pot_rhc = 0;
};

double nice_up(double x, double step)
{
    if (x <= 0) return step;
    return std::ceil(x / step) * step;
}

// "Nice" log-scale bounds using 1/2/5 decades.
// These are chosen to bracket the data with a bit of padding.
double nice_log_down(double x)
{
    if (x <= 0.0) return 0.1;
    const double e = std::floor(std::log10(x));
    const double base = std::pow(10.0, e);
    const double m = x / base; // in [1,10)

    double nice_m = 1.0;
    if (m >= 5.0) nice_m = 5.0;
    else if (m >= 2.0) nice_m = 2.0;
    else nice_m = 1.0;

    return nice_m * base;
}

double nice_log_up(double x)
{
    if (x <= 0.0) return 1.0;
    const double e = std::floor(std::log10(x));
    const double base = std::pow(10.0, e);
    const double m = x / base; // in [1,10)

    double nice_m = 1.0;
    if (m <= 1.0) nice_m = 1.0;
    else if (m <= 2.0) nice_m = 2.0;
    else if (m <= 5.0) nice_m = 5.0;
    else nice_m = 10.0;

    return nice_m * base;
}

// Map a positive power value (kW) onto left-axis coordinates [0, y_pot_max]
// using a log10 scale defined by [pmin, pmax].
double power_kw_to_left_axis_log(double pkw, double pmin, double pmax, double y_pot_max)
{
    if (!(pmin > 0.0) || !(pmax > pmin) || !(y_pot_max > 0.0)) return 0.0;
    if (pkw <= pmin) return 0.0;
    if (pkw >= pmax) return y_pot_max;

    const double lmin = std::log10(pmin);
    const double lmax = std::log10(pmax);
    const double lval = std::log10(pkw);
    const double frac = (lval - lmin) / (lmax - lmin);
    return frac * y_pot_max;
}

plot_data compute_plot_data(const TH1D &h_bnb, const TH1D &h_fhc, const TH1D &h_rhc)
{
    plot_data d;
    const int nbins = h_bnb.GetNbinsX();

    // Step arrays: (xlo,0) then for each bin: (x_right, y_prev) and (x_right, y_new)
    d.x_step.resize(1 + 2 * nbins);
    d.y_total.resize(1 + 2 * nbins);
    d.y_bnb.resize(1 + 2 * nbins);
    d.y_fhc.resize(1 + 2 * nbins);
    d.y_rhc.resize(1 + 2 * nbins);

    const double xlo = h_bnb.GetXaxis()->GetXmin();
    d.x_step[0] = xlo;
    d.y_total[0] = 0;
    d.y_bnb[0] = 0;
    d.y_fhc[0] = 0;
    d.y_rhc[0] = 0;

    // Energies (GeV) for power computation.
    // BNB: 8 GeV Booster
    // NuMI: 120 GeV Main Injector
    const double E_BNB_GeV = 8.0;
    const double E_NUMI_GeV = 120.0;

    // Conversion factor: GeV -> J, then /1000 for kW.
    const double GeV_to_kW_s = 1.602176634e-13; // (J/GeV)/1000

    double cum_total = 0;
    double cum_bnb = 0;
    double cum_fhc = 0;
    double cum_rhc = 0;

    double max_cum = 0;
    double max_pow = 0;
    double min_pow = std::numeric_limits<double>::infinity();

    int idx = 1;
    for (int i = 1; i <= nbins; ++i)
    {
        const double pot_bnb = h_bnb.GetBinContent(i);
        const double pot_fhc = h_fhc.GetBinContent(i);
        const double pot_rhc = h_rhc.GetBinContent(i);
        const double pot_tot = pot_bnb + pot_fhc + pot_rhc;

        d.pot_bnb += pot_bnb;
        d.pot_fhc += pot_fhc;
        d.pot_rhc += pot_rhc;
        d.pot_total += pot_tot;

        const double x_right = h_bnb.GetXaxis()->GetBinUpEdge(i);
        const double x_center = h_bnb.GetXaxis()->GetBinCenter(i);
        const double binw = h_bnb.GetXaxis()->GetBinWidth(i); // seconds (since x is unix seconds)

        // Horizontal to the right edge at previous cumulative.
        d.x_step[idx] = x_right;
        d.y_total[idx] = cum_total;
        d.y_bnb[idx] = cum_bnb;
        d.y_fhc[idx] = cum_fhc;
        d.y_rhc[idx] = cum_rhc;
        ++idx;

        // Update cumulative (in units of 1e20 POT).
        cum_bnb += pot_bnb / 1e20;
        cum_fhc += pot_fhc / 1e20;
        cum_rhc += pot_rhc / 1e20;
        cum_total += pot_tot / 1e20;

        // Vertical at the right edge to the new cumulative.
        d.x_step[idx] = x_right;
        d.y_total[idx] = cum_total;
        d.y_bnb[idx] = cum_bnb;
        d.y_fhc[idx] = cum_fhc;
        d.y_rhc[idx] = cum_rhc;
        ++idx;

        max_cum = std::max(max_cum, cum_total);

        // Beam power points (kW), averaged over the bin width.
        if (binw > 0)
        {
            if (pot_bnb > 0)
            {
                const double pkw = (pot_bnb / binw) * (E_BNB_GeV * GeV_to_kW_s);
                d.x_pow_bnb.push_back(x_center);
                d.y_pow_bnb_scaled.push_back(pkw);
                max_pow = std::max(max_pow, pkw);
                min_pow = std::min(min_pow, pkw);
            }
            if (pot_fhc > 0)
            {
                const double pkw = (pot_fhc / binw) * (E_NUMI_GeV * GeV_to_kW_s);
                d.x_pow_fhc.push_back(x_center);
                d.y_pow_fhc_scaled.push_back(pkw);
                max_pow = std::max(max_pow, pkw);
                min_pow = std::min(min_pow, pkw);
            }
            if (pot_rhc > 0)
            {
                const double pkw = (pot_rhc / binw) * (E_NUMI_GeV * GeV_to_kW_s);
                d.x_pow_rhc.push_back(x_center);
                d.y_pow_rhc_scaled.push_back(pkw);
                max_pow = std::max(max_pow, pkw);
                min_pow = std::min(min_pow, pkw);
            }
        }
    }

    // Nice left-axis max (like the reference: 0,5,10,...)
    d.y_pot_max = nice_up(max_cum, 5.0) * 1.05;

    // Log-scale right-axis bounds (kW). Add a little padding so points do not sit on the frame.
    if (std::isfinite(min_pow) && max_pow > 0.0)
    {
        const double pad_lo = min_pow / 1.25;
        const double pad_hi = max_pow * 1.25;
        d.pow_axis_min = nice_log_down(pad_lo);
        d.pow_axis_max = nice_log_up(pad_hi);
        if (!(d.pow_axis_min > 0.0)) d.pow_axis_min = nice_log_down(min_pow);
        if (!(d.pow_axis_max > d.pow_axis_min)) d.pow_axis_max = nice_log_up(max_pow);
        if (!(d.pow_axis_max > d.pow_axis_min))
        {
            // Fallback (should not happen for realistic inputs)
            d.pow_axis_min = 1.0;
            d.pow_axis_max = 1000.0;
        }
    }
    else
    {
        // No power points: pick a harmless default log range.
        d.pow_axis_min = 1.0;
        d.pow_axis_max = 1000.0;
    }

    // Map power arrays onto left-axis coordinates using log10 scaling.
    for (double &pkw : d.y_pow_bnb_scaled)
    {
        pkw = power_kw_to_left_axis_log(pkw, d.pow_axis_min, d.pow_axis_max, d.y_pot_max);
    }
    for (double &pkw : d.y_pow_fhc_scaled)
    {
        pkw = power_kw_to_left_axis_log(pkw, d.pow_axis_min, d.pow_axis_max, d.y_pot_max);
    }
    for (double &pkw : d.y_pow_rhc_scaled)
    {
        pkw = power_kw_to_left_axis_log(pkw, d.pow_axis_min, d.pow_axis_max, d.y_pot_max);
    }

    return d;
}

// -------------------- Drawing --------------------

void draw_plot(const TH1D &h_bnb, const TH1D &h_fhc, const TH1D &h_rhc,
               const plot_data &d,
               const std::vector<run_period> &runs,
               const time_t tmin, const time_t tmax,
               const std::string &outfile)
{
    TCanvas canvas("c", "POT timeline",
                   kCanvasWidth,
                   kCanvasHeight);

    // Two pads: top legend strip + main plot.
    const double ml = 0.08;
    const double mr = 0.14; // extra room for log-scale RHS labels (e.g. 10^{3})
    const double split = 0.82;
    TPad *main_pad = new TPad("pad_main", "pad_main", 0.0, 0.0, 1.0, split);
    TPad *legend_pad = new TPad("pad_legend", "pad_legend", 0.0, split, 1.0, 1.0);

    main_pad->SetTopMargin(0.02);
    main_pad->SetBottomMargin(0.20);
    main_pad->SetLeftMargin(ml);
    main_pad->SetRightMargin(mr);

    legend_pad->SetTopMargin(0.10);
    legend_pad->SetBottomMargin(0.02);
    legend_pad->SetLeftMargin(ml);
    legend_pad->SetRightMargin(mr);

    main_pad->Draw();
    legend_pad->Draw();

    // ---- Main plot ----
    main_pad->cd();
    main_pad->SetTickx(1);
    main_pad->SetTicky(0);
    main_pad->SetGridx(true);
    main_pad->SetGridy(true);

    const double xlo = h_bnb.GetXaxis()->GetXmin();
    const double xhi = h_bnb.GetXaxis()->GetXmax();

    // Frame
    TH1F *frame = main_pad->DrawFrame(xlo, 0.0, xhi, d.y_pot_max);
    frame->GetXaxis()->SetTimeDisplay(1);
    frame->GetXaxis()->SetTimeOffset(0, "gmt");
    frame->GetXaxis()->SetTimeFormat("%Y");
    // Avoid duplicate year labels like "2016 2016" when ROOT optimizes time ticks.
    // Keep fewer major divisions for a year-only axis.
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitle("Year");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleOffset(1.10);

    frame->GetYaxis()->SetTitle("Accumulated POT (#times 10^{20})");
    frame->GetYaxis()->SetNdivisions(507);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleOffset(0.70);

    // Run bands (allocate on heap so they survive through SaveAs)
    const Int_t band_col = TColor::GetColor("#f8bbd0"); // light pink
    std::vector<TBox *> boxes;
    boxes.reserve(runs.size());

    for (size_t i = 0; i < runs.size(); ++i)
    {
        const float alpha = (i % 2 == 0) ? 0.14f : 0.08f;
        TBox *b = new TBox(runs[i].x1, 0.0, runs[i].x2, d.y_pot_max);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 22, 0)
        b->SetFillColorAlpha(band_col, alpha);
#else
        b->SetFillColor(band_col);
#endif
        b->SetLineColor(0);
        b->Draw("same");
        boxes.push_back(b);

        // Label near top
        TLatex tl;
        tl.SetTextFont(42);
        tl.SetTextSize(0.040);
        tl.SetTextAlign(22); // center
        const double xm = 0.5 * (runs[i].x1 + runs[i].x2);
        const double ym = 0.975 * d.y_pot_max;
        tl.DrawLatex(xm, ym, runs[i].label.c_str());
    }

    // Colors
    const Int_t col_total = TColor::GetColor("#1e88e5"); // blue
    const Int_t col_bnb   = TColor::GetColor("#00c853"); // green
    const Int_t col_fhc   = TColor::GetColor("#ffb300"); // amber
    const Int_t col_rhc   = TColor::GetColor("#d81b60"); // magenta-ish

    // Solid lines for all cumulative curves (use colour + underlay for separation).
    const Int_t ls_total = 1; // solid
    const Int_t ls_bnb   = 1; // solid
    const Int_t ls_fhc   = 1; // solid
    const Int_t ls_rhc   = 1; // solid

    // Cumulative step curves
    TGraph g_total(d.x_step.size(), d.x_step.data(), d.y_total.data());
    TGraph g_bnb(d.x_step.size(), d.x_step.data(), d.y_bnb.data());
    TGraph g_fhc(d.x_step.size(), d.x_step.data(), d.y_fhc.data());
    TGraph g_rhc(d.x_step.size(), d.x_step.data(), d.y_rhc.data());

    // Outline for contrast (black underlay), like your existing macro
    TGraph g_total_outline(d.x_step.size(), d.x_step.data(), d.y_total.data());
    TGraph g_bnb_outline(d.x_step.size(), d.x_step.data(), d.y_bnb.data());
    TGraph g_fhc_outline(d.x_step.size(), d.x_step.data(), d.y_fhc.data());
    TGraph g_rhc_outline(d.x_step.size(), d.x_step.data(), d.y_rhc.data());

    g_total_outline.SetLineColor(kBlack);
    g_total_outline.SetLineStyle(ls_total);
    g_total_outline.SetLineWidth(3);
    g_bnb_outline.SetLineColor(kBlack);
    g_bnb_outline.SetLineStyle(ls_bnb);
    g_bnb_outline.SetLineWidth(3);
    g_fhc_outline.SetLineColor(kBlack);
    g_fhc_outline.SetLineStyle(ls_fhc);
    g_fhc_outline.SetLineWidth(3);
    g_rhc_outline.SetLineColor(kBlack);
    g_rhc_outline.SetLineStyle(ls_rhc);
    g_rhc_outline.SetLineWidth(3);

    g_total.SetLineColor(col_total);
    g_total.SetLineStyle(ls_total);
    g_total.SetLineWidth(2);
    g_bnb.SetLineColor(col_bnb);
    g_bnb.SetLineStyle(ls_bnb);
    g_bnb.SetLineWidth(2);
    g_fhc.SetLineColor(col_fhc);
    g_fhc.SetLineStyle(ls_fhc);
    g_fhc.SetLineWidth(2);
    g_rhc.SetLineColor(col_rhc);
    g_rhc.SetLineStyle(ls_rhc);
    g_rhc.SetLineWidth(2);

    g_total_outline.Draw("L SAME");
    g_bnb_outline.Draw("L SAME");
    g_fhc_outline.Draw("L SAME");
    g_rhc_outline.Draw("L SAME");

    g_total.Draw("L SAME");
    g_bnb.Draw("L SAME");
    g_fhc.Draw("L SAME");
    g_rhc.Draw("L SAME");

    // Beam power points (already mapped onto left axis)
    TGraph gp_bnb(d.x_pow_bnb.size(), d.x_pow_bnb.data(), d.y_pow_bnb_scaled.data());
    TGraph gp_fhc(d.x_pow_fhc.size(), d.x_pow_fhc.data(), d.y_pow_fhc_scaled.data());
    TGraph gp_rhc(d.x_pow_rhc.size(), d.x_pow_rhc.data(), d.y_pow_rhc_scaled.data());

    gp_bnb.SetLineColor(col_bnb);
    gp_bnb.SetMarkerColor(col_bnb);
    gp_bnb.SetMarkerStyle(24); // open markers read better over run bands
    gp_bnb.SetMarkerSize(0.65);
    gp_fhc.SetLineColor(col_fhc);
    gp_fhc.SetMarkerColor(col_fhc);
    gp_fhc.SetMarkerStyle(25);
    gp_fhc.SetMarkerSize(0.65);
    gp_rhc.SetLineColor(col_rhc);
    gp_rhc.SetMarkerColor(col_rhc);
    gp_rhc.SetMarkerStyle(26);
    gp_rhc.SetMarkerSize(0.65);

    gp_bnb.Draw("P SAME");
    gp_fhc.Draw("P SAME");
    gp_rhc.Draw("P SAME");

    // Right axis: beam power (kW), logarithmic.
    //  - "+"  : tick marks on positive side (into the right margin)
    //  - "L"  : left-adjust labels so they extend rightwards (outside the frame)
    //  - "G"  : loGarithmic scale
    //  - "S"  : tick size uses fTickSize (SetTickSize)
    TGaxis right_axis(xhi, 0.0, xhi, d.y_pot_max, d.pow_axis_min, d.pow_axis_max, 510, "+LGS");
    right_axis.SetLineColor(kRed + 1);
    right_axis.SetLabelColor(kRed + 1);
    right_axis.SetTitleColor(kRed + 1);
    right_axis.SetLineWidth(2);
    right_axis.SetLabelFont(42);
    right_axis.SetTitleFont(42);
    right_axis.SetLabelSize(0.045);
    right_axis.SetTitleSize(0.045);
    right_axis.SetTitleOffset(0.85);
    right_axis.SetTickSize(0.02);
    right_axis.SetTitle("Beam Power (kW)");
    right_axis.Draw();

    // Small date range + totals, in the bottom margin (NDC in pad coordinates)
    {
        const double total = d.pot_total;
        const double bnb   = d.pot_bnb;
        const double fhc   = d.pot_fhc;
        const double rhc   = d.pot_rhc;

        TLatex infoL;
        infoL.SetNDC(true);
        infoL.SetTextFont(42);
        infoL.SetTextSize(0.038);
        infoL.SetTextAlign(13); // left, top

        const std::string s0 = utc_date_ymd(tmin);
        const std::string s1 = utc_date_ymd(tmax);

        const double y0 = 0.025;
        infoL.DrawLatex(ml + 0.01, y0 + 0.080, Form("%s  --  %s", s0.c_str(), s1.c_str()));
        infoL.DrawLatex(ml + 0.01, y0 + 0.030, Form("Total POT: %.3g", total));
        if (total > 0)
        {
            // Right-anchored so the NuMI breakdown does not get clipped by the pad edge.
            TLatex infoR;
            infoR.SetNDC(true);
            infoR.SetTextFont(42);
            infoR.SetTextSize(0.038);
            infoR.SetTextAlign(33); // right, top
            const double xr = 1.0 - mr - 0.01;
            infoR.DrawLatex(xr, y0 + 0.080, Form("BNB: %.3g (%.1f%%)", bnb, 100.0 * bnb / total));
            infoR.DrawLatex(xr, y0 + 0.030, Form("NuMI FHC: %.3g (%.1f%%)   NuMI RHC: %.3g (%.1f%%)",
                                                 fhc, 100.0 * fhc / total, rhc, 100.0 * rhc / total));
        }
    }

    main_pad->RedrawAxis();

    // ---- Legend strip ----
    legend_pad->cd();
    const double lx1 = legend_pad->GetLeftMargin() + 0.02;
    const double lx2 = 1.0 - legend_pad->GetRightMargin() - 0.02;
    const double ly1 = legend_pad->GetBottomMargin() + 0.02;
    const double ly2 = 1.0 - legend_pad->GetTopMargin() - 0.02;

    TLegend legend(lx1, ly1, lx2, ly2);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetNColumns(2);
    legend.SetColumnSeparation(0.10);
    legend.SetEntrySeparation(0.02);
    legend.SetMargin(0.22);

    const double s_main = 0.045;
    const double s_leg = s_main * (split / (1.0 - split));
    legend.SetTextSize(s_leg);

    legend.AddEntry(&g_total, "Total accumulated POT", "l");
    legend.AddEntry(&g_bnb,   "BNB accumulated POT", "l");
    legend.AddEntry(&g_fhc,   "NuMI-FHC accumulated POT", "l");
    legend.AddEntry(&g_rhc,   "NuMI-RHC accumulated POT", "l");

    legend.AddEntry(&gp_bnb, "BNB beam power (avg/bin)", "p");
    legend.AddEntry(&gp_fhc, "NuMI-FHC beam power (avg/bin)", "p");
    legend.AddEntry(&gp_rhc, "NuMI-RHC beam power (avg/bin)", "p");
    legend.Draw();

    // Save
    canvas.SaveAs(outfile.c_str());

    // Cleanup (pads). Boxes are owned by ROOT pads after Draw(); OK to leak in a macro,
    // but delete the pads (as in your existing macro).
    delete main_pad;
    delete legend_pad;
}

// -------------------- Main macro entry --------------------

void plotPotTimelineStyleInternal(const char *out = nullptr)
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

    // Same queries as your plotPotSimple.C
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
    time_t tmin = 0;
    time_t tmax = 0;

    const std::vector<pot_samples> all_samples{bnb_samples, fhc_samples, rhc_samples};
    if (!determine_range_weekly(all_samples, xlo, xhi, nbins, tmin, tmax))
    {
        std::cerr << "No time range\n";
        return;
    }

    // Weekly POT histograms (in POT, not scaled)
    TH1D hBNB("hBNB", "", nbins, xlo, xhi);
    TH1D hFHC("hFHC", "", nbins, xlo, xhi);
    TH1D hRHC("hRHC", "", nbins, xlo, xhi);
    hBNB.SetDirectory(nullptr);
    hFHC.SetDirectory(nullptr);
    hRHC.SetDirectory(nullptr);

    fill_histogram_pot(bnb_samples, hBNB);
    fill_histogram_pot(fhc_samples, hFHC);
    fill_histogram_pot(rhc_samples, hRHC);

    // Build plot data
    const plot_data d = compute_plot_data(hBNB, hFHC, hRHC);

    // Run bands inferred from periods with significant total POT in a week.
    // Tune these if you want fewer/more run blocks.
    const bool draw_runs = true;
    const double active_threshold_pot = 1e17;           // "active" if > 1e17 POT in the bin
    const double merge_gap_seconds = 42.0 * 86400.0;    // merge gaps shorter than 42 days

    std::vector<run_period> runs;
    if (draw_runs)
    {
        runs = infer_run_periods(hBNB, hFHC, hRHC, active_threshold_pot, merge_gap_seconds);
    }

    const auto out_path = resolve_output_file(out ? out : "", "pot_timeline_style");
    const std::string outfile = out_path.string();

    draw_plot(hBNB, hFHC, hRHC, d, runs, tmin, tmax, outfile);
}


int nuxsec_plot(const char *out = nullptr)
{
    plotPotTimelineStyleInternal(out);
    return 0;
}

void plotPotTimelineStyle(const char *out = nullptr)
{
    plotPotTimelineStyleInternal(out);
}
