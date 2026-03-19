// plot_score_systematic_envelopes.C
//
// Build first-score systematic-envelope diagnostics using the compressed universe branches
// written by EventWeightAnalysis.
//
// For each available family among
//   - weightsGenie
//   - weightsFlux
//   - weightsReint
//   - weightsPPFX
//
// the macro builds, after a base selection:
//
//   1) the CV score histogram for truth signal,
//   2) the CV score histogram for total background,
//   3) the corresponding universe-weighted histograms,
//   4) the cumulative tails N(score > t),
//   5) ratio-to-CV bands per score bin and per tail threshold.
//
// Scores below xmin are dropped; scores above xmax are accumulated in the last bin so
// the cumulative tail keeps the overflow contribution.
//
// By default the total-background envelope varies only the MC component; EXT is kept fixed in
// every universe and added to both numerator and denominator of the background ratio.
//
// Notes on weights:
//   - The universe branches are decoded from unsigned short as w = raw / 1000.
//   - The varied weight is taken as
//         w_var = w_base * w_universe / w_family_cv
//     where w_family_cv defaults to 1.0, except PPFX where the default expression is "ppfx_cv".
//   - If your base weight does NOT already contain ppfx_cv, set ppfx_cv_expr="" when calling.
//
// Example:
//   ./heron macro macros/plot_score_systematic_envelopes.C
//
//   ./heron macro macros/plot_score_systematic_envelopes.C \
//     'plot_score_systematic_envelopes("./scratch/out/event_list_myana.root",
//                                     "sel_muon",
//                                     "is_signal",
//                                     "w_nominal",
//                                     50, -15, 15,
//                                     true,
//                                     "both",
//                                     "ppfx_cv",
//                                     "first_score_syst_env")'

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "PlottingHelper.hh"
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

bool is_simple_identifier(const std::string &expr)
{
    if (expr.empty())
        return false;

    const unsigned char first = static_cast<unsigned char>(expr.front());
    if (!(std::isalpha(first) || expr.front() == '_'))
        return false;

    for (std::size_t i = 1; i < expr.size(); ++i)
    {
        const unsigned char c = static_cast<unsigned char>(expr[i]);
        if (!(std::isalnum(c) || expr[i] == '_'))
            return false;
    }

    return true;
}

ROOT::RDF::RNode filter_by_mask(ROOT::RDF::RNode node,
                                std::shared_ptr<const std::vector<char>> mask,
                                const std::string &sample_id_col = "sample_id")
{
    return node.Filter(
        [mask](int sid) {
            return sid >= 0 &&
                   sid < static_cast<int>(mask->size()) &&
                   (*mask)[static_cast<std::size_t>(sid)];
        },
        {sample_id_col});
}

double decode_packed_weight(unsigned short raw)
{
    return 1.0e-3 * static_cast<double>(raw);
}

std::vector<double> reverse_cumulative(const std::vector<double> &bins)
{
    std::vector<double> out(bins.size(), 0.0);
    double running = 0.0;
    for (int i = static_cast<int>(bins.size()) - 1; i >= 0; --i)
    {
        running += bins[static_cast<std::size_t>(i)];
        out[static_cast<std::size_t>(i)] = running;
    }
    return out;
}

std::vector<std::vector<double>> reverse_cumulative_set(const std::vector<std::vector<double>> &hset)
{
    std::vector<std::vector<double>> out;
    out.reserve(hset.size());
    for (const auto &h : hset)
        out.push_back(reverse_cumulative(h));
    return out;
}

struct FamilyConfig
{
    std::string key;
    std::string title;
    std::string branch;
    std::string cv_expr;
};

struct EnvelopeAccumulator
{
    int nbins = 0;
    double xmin = 0.0;
    double xmax = 1.0;

    std::vector<double> sig_cv;
    std::vector<double> bkg_mc_cv;
    std::vector<double> bkg_ext_cv;

    std::vector<std::vector<double>> sig_univ;
    std::vector<std::vector<double>> bkg_mc_univ;

    std::size_t nuniv = 0;
    bool warned_size_mismatch = false;

    EnvelopeAccumulator() = default;

    EnvelopeAccumulator(int nbin, double xlo, double xhi)
        : nbins(nbin), xmin(xlo), xmax(xhi),
          sig_cv(static_cast<std::size_t>(nbin), 0.0),
          bkg_mc_cv(static_cast<std::size_t>(nbin), 0.0),
          bkg_ext_cv(static_cast<std::size_t>(nbin), 0.0)
    {
    }

    int find_bin(double x) const
    {
        if (!std::isfinite(x))
            return -1;
        if (x < xmin)
            return -1;
        if (x >= xmax)
            return nbins - 1;

        const double frac = (x - xmin) / (xmax - xmin);
        int ibin = static_cast<int>(std::floor(frac * static_cast<double>(nbins)));
        if (ibin < 0)
            ibin = 0;
        if (ibin >= nbins)
            ibin = nbins - 1;
        return ibin;
    }

    void initialise_universes(std::size_t n)
    {
        nuniv = n;
        sig_univ.assign(n, std::vector<double>(static_cast<std::size_t>(nbins), 0.0));
        bkg_mc_univ.assign(n, std::vector<double>(static_cast<std::size_t>(nbins), 0.0));
    }

    void ensure_universe_size(std::size_t n, const std::string &family_key)
    {
        if (n == 0)
            return;

        if (nuniv == 0)
        {
            initialise_universes(n);
            return;
        }

        if (n != nuniv && !warned_size_mismatch)
        {
            std::cerr << "[plot_score_systematic_envelopes] warning: family '" << family_key
                      << "' saw universe-vector size mismatch. First size=" << nuniv
                      << ", later size=" << n
                      << ". Using the first size and treating missing entries as unity.\n";
            warned_size_mismatch = true;
        }
    }

    void fill_mc(bool is_signal,
                 double score,
                 double w_base,
                 double family_cv,
                 const ROOT::RVec<unsigned short> &packed_universes,
                 const std::string &family_key)
    {
        const int ibin = find_bin(score);
        if (ibin < 0)
            return;
        if (!std::isfinite(w_base))
            return;

        if (nuniv == 0 && !packed_universes.empty())
        {
            initialise_universes(packed_universes.size());
            for (std::size_t iu = 0; iu < nuniv; ++iu)
            {
                sig_univ[iu] = sig_cv;
                bkg_mc_univ[iu] = bkg_mc_cv;
            }
        }
        else
        {
            ensure_universe_size(packed_universes.size(), family_key);
        }

        if (is_signal)
            sig_cv[static_cast<std::size_t>(ibin)] += w_base;
        else
            bkg_mc_cv[static_cast<std::size_t>(ibin)] += w_base;

        if (nuniv == 0)
            return;

        double denom = 1.0;
        if (std::isfinite(family_cv) && family_cv > 0.0)
            denom = family_cv;

        for (std::size_t iu = 0; iu < nuniv; ++iu)
        {
            double w_var = w_base;
            if (iu < packed_universes.size())
            {
                double w_family = decode_packed_weight(packed_universes[iu]);
                if (!std::isfinite(w_family) || w_family < 0.0)
                    w_family = 0.0;
                w_var = w_base * w_family / denom;
            }

            if (!std::isfinite(w_var))
                continue;

            if (is_signal)
                sig_univ[iu][static_cast<std::size_t>(ibin)] += w_var;
            else
                bkg_mc_univ[iu][static_cast<std::size_t>(ibin)] += w_var;
        }
    }

    void fill_ext(double score, double w_base)
    {
        const int ibin = find_bin(score);
        if (ibin < 0)
            return;
        if (!std::isfinite(w_base))
            return;
        bkg_ext_cv[static_cast<std::size_t>(ibin)] += w_base;
    }

    void merge(const EnvelopeAccumulator &other)
    {
        if (nbins != other.nbins || xmin != other.xmin || xmax != other.xmax)
            return;

        if (nuniv == 0 && other.nuniv > 0)
            initialise_universes(other.nuniv);

        if (other.nuniv != 0 && nuniv != other.nuniv && !warned_size_mismatch)
        {
            std::cerr << "[plot_score_systematic_envelopes] warning: slot-merge universe-size mismatch. "
                      << "Keeping first size=" << nuniv << ", other size=" << other.nuniv << ".\n";
            warned_size_mismatch = true;
        }

        for (int b = 0; b < nbins; ++b)
        {
            const std::size_t ib = static_cast<std::size_t>(b);
            sig_cv[ib] += other.sig_cv[ib];
            bkg_mc_cv[ib] += other.bkg_mc_cv[ib];
            bkg_ext_cv[ib] += other.bkg_ext_cv[ib];
        }

        const std::size_t nmerge = std::min(nuniv, other.nuniv);
        for (std::size_t iu = 0; iu < nmerge; ++iu)
        {
            for (int b = 0; b < nbins; ++b)
            {
                const std::size_t ib = static_cast<std::size_t>(b);
                sig_univ[iu][ib] += other.sig_univ[iu][ib];
                bkg_mc_univ[iu][ib] += other.bkg_mc_univ[iu][ib];
            }
        }
    }

    std::vector<double> total_bkg_cv(bool include_ext) const
    {
        std::vector<double> out = bkg_mc_cv;
        if (include_ext)
        {
            for (std::size_t i = 0; i < out.size(); ++i)
                out[i] += bkg_ext_cv[i];
        }
        return out;
    }

    std::vector<std::vector<double>> total_bkg_univ(bool include_ext) const
    {
        std::vector<std::vector<double>> out = bkg_mc_univ;
        if (include_ext)
        {
            for (auto &u : out)
            {
                for (std::size_t i = 0; i < u.size(); ++i)
                    u[i] += bkg_ext_cv[i];
            }
        }
        return out;
    }
};

struct EnvelopeSummary
{
    int nbins = 0;
    double xmin = 0.0;
    double xmax = 1.0;
    std::size_t nuniv = 0;

    std::vector<double> cv;
    std::vector<double> mean_ratio;
    std::vector<double> rms_ratio;
    std::vector<double> min_ratio;
    std::vector<double> max_ratio;

    double cv_max = 0.0;
    double cv_min_positive = std::numeric_limits<double>::infinity();
    double ratio_min_global = 0.95;
    double ratio_max_global = 1.05;
};

EnvelopeSummary summarise_envelope(const std::vector<double> &cv,
                                   const std::vector<std::vector<double>> &universes,
                                   double xmin,
                                   double xmax)
{
    EnvelopeSummary out;
    out.nbins = static_cast<int>(cv.size());
    out.xmin = xmin;
    out.xmax = xmax;
    out.nuniv = universes.size();
    out.cv = cv;
    out.mean_ratio.assign(cv.size(), std::numeric_limits<double>::quiet_NaN());
    out.rms_ratio.assign(cv.size(), std::numeric_limits<double>::quiet_NaN());
    out.min_ratio.assign(cv.size(), std::numeric_limits<double>::quiet_NaN());
    out.max_ratio.assign(cv.size(), std::numeric_limits<double>::quiet_NaN());

    bool have_ratio = false;

    for (std::size_t ib = 0; ib < cv.size(); ++ib)
    {
        const double cv_bin = cv[ib];
        out.cv_max = std::max(out.cv_max, cv_bin);
        if (cv_bin > 0.0)
            out.cv_min_positive = std::min(out.cv_min_positive, cv_bin);

        if (!(cv_bin > 0.0) || universes.empty())
            continue;

        double sum = 0.0;
        double sumsq_about_cv = 0.0;
        double rmin = std::numeric_limits<double>::infinity();
        double rmax = -std::numeric_limits<double>::infinity();
        int nvalid = 0;

        for (const auto &u : universes)
        {
            if (ib >= u.size())
                continue;

            const double r = u[ib] / cv_bin;
            if (!std::isfinite(r))
                continue;

            sum += r;
            sumsq_about_cv += (r - 1.0) * (r - 1.0);
            rmin = std::min(rmin, r);
            rmax = std::max(rmax, r);
            ++nvalid;
        }

        if (nvalid == 0)
            continue;

        const double mean = sum / static_cast<double>(nvalid);
        const double rms = std::sqrt(sumsq_about_cv / static_cast<double>(nvalid));

        out.mean_ratio[ib] = mean;
        out.rms_ratio[ib] = rms;
        out.min_ratio[ib] = rmin;
        out.max_ratio[ib] = rmax;

        out.ratio_min_global = std::min(out.ratio_min_global, std::min(rmin, 1.0 - rms));
        out.ratio_max_global = std::max(out.ratio_max_global, std::max(rmax, 1.0 + rms));
        have_ratio = true;
    }

    if (!have_ratio)
    {
        out.ratio_min_global = 0.95;
        out.ratio_max_global = 1.05;
    }

    return out;
}

std::unique_ptr<TH1D> make_histogram(const std::string &name,
                                     int nbins,
                                     double xmin,
                                     double xmax,
                                     const std::vector<double> &contents)
{
    auto h = std::make_unique<TH1D>(name.c_str(), "", nbins, xmin, xmax);
    h->SetDirectory(nullptr);
    for (int b = 1; b <= nbins; ++b)
    {
        const std::size_t ib = static_cast<std::size_t>(b - 1);
        if (ib < contents.size())
            h->SetBinContent(b, contents[ib]);
    }
    return h;
}

std::unique_ptr<TGraphAsymmErrors> make_rms_band(const EnvelopeSummary &s, int color)
{
    const double bw = (s.xmax - s.xmin) / static_cast<double>(s.nbins);

    std::vector<double> x, y, exl, exh, eyl, eyh;
    x.reserve(static_cast<std::size_t>(s.nbins));
    y.reserve(static_cast<std::size_t>(s.nbins));
    exl.reserve(static_cast<std::size_t>(s.nbins));
    exh.reserve(static_cast<std::size_t>(s.nbins));
    eyl.reserve(static_cast<std::size_t>(s.nbins));
    eyh.reserve(static_cast<std::size_t>(s.nbins));

    for (int b = 0; b < s.nbins; ++b)
    {
        const std::size_t ib = static_cast<std::size_t>(b);
        if (!(s.cv[ib] > 0.0) || !std::isfinite(s.rms_ratio[ib]))
            continue;

        x.push_back(s.xmin + (static_cast<double>(b) + 0.5) * bw);
        y.push_back(1.0);
        exl.push_back(0.5 * bw);
        exh.push_back(0.5 * bw);
        eyl.push_back(s.rms_ratio[ib]);
        eyh.push_back(s.rms_ratio[ib]);
    }

    auto g = std::make_unique<TGraphAsymmErrors>(static_cast<int>(x.size()),
                                                 x.data(), y.data(),
                                                 exl.data(), exh.data(),
                                                 eyl.data(), eyh.data());
    g->SetFillColorAlpha(color, 0.28);
    g->SetLineColorAlpha(color, 0.0);
    g->SetMarkerSize(0.0);
    return g;
}

std::unique_ptr<TGraphAsymmErrors> make_env_band(const EnvelopeSummary &s, int color)
{
    const double bw = (s.xmax - s.xmin) / static_cast<double>(s.nbins);

    std::vector<double> x, y, exl, exh, eyl, eyh;
    x.reserve(static_cast<std::size_t>(s.nbins));
    y.reserve(static_cast<std::size_t>(s.nbins));
    exl.reserve(static_cast<std::size_t>(s.nbins));
    exh.reserve(static_cast<std::size_t>(s.nbins));
    eyl.reserve(static_cast<std::size_t>(s.nbins));
    eyh.reserve(static_cast<std::size_t>(s.nbins));

    for (int b = 0; b < s.nbins; ++b)
    {
        const std::size_t ib = static_cast<std::size_t>(b);
        if (!(s.cv[ib] > 0.0) || !std::isfinite(s.min_ratio[ib]) || !std::isfinite(s.max_ratio[ib]))
            continue;

        const double ymid = 0.5 * (s.min_ratio[ib] + s.max_ratio[ib]);
        x.push_back(s.xmin + (static_cast<double>(b) + 0.5) * bw);
        y.push_back(ymid);
        exl.push_back(0.5 * bw);
        exh.push_back(0.5 * bw);
        eyl.push_back(ymid - s.min_ratio[ib]);
        eyh.push_back(s.max_ratio[ib] - ymid);
    }

    auto g = std::make_unique<TGraphAsymmErrors>(static_cast<int>(x.size()),
                                                 x.data(), y.data(),
                                                 exl.data(), exh.data(),
                                                 eyl.data(), eyh.data());
    g->SetFillColorAlpha(color, 0.18);
    g->SetLineColorAlpha(color, 0.0);
    g->SetMarkerSize(0.0);
    return g;
}

std::unique_ptr<TGraph> make_ratio_curve(const EnvelopeSummary &s,
                                         const std::vector<double> &yvals,
                                         int color,
                                         int linestyle,
                                         int linewidth)
{
    const double bw = (s.xmax - s.xmin) / static_cast<double>(s.nbins);

    std::vector<double> x, y;
    x.reserve(static_cast<std::size_t>(s.nbins));
    y.reserve(static_cast<std::size_t>(s.nbins));

    for (int b = 0; b < s.nbins; ++b)
    {
        const std::size_t ib = static_cast<std::size_t>(b);
        if (!(s.cv[ib] > 0.0) || !std::isfinite(yvals[ib]))
            continue;

        x.push_back(s.xmin + (static_cast<double>(b) + 0.5) * bw);
        y.push_back(yvals[ib]);
    }

    auto g = std::make_unique<TGraph>(static_cast<int>(x.size()), x.data(), y.data());
    g->SetLineColor(color);
    g->SetLineStyle(linestyle);
    g->SetLineWidth(linewidth);
    g->SetMarkerSize(0.0);
    return g;
}

double sum_contents(const std::vector<double> &v)
{
    double out = 0.0;
    for (double x : v)
        out += x;
    return out;
}

void draw_single_envelope_plot(const EnvelopeSummary &s,
                               const std::string &canvas_name,
                               const std::string &title_text,
                               const std::string &x_title,
                               const std::string &y_title_top,
                               const std::string &output_stem,
                               int line_color,
                               bool logy_top,
                               const std::string &envelope_mode)
{
    if (s.nuniv == 0)
    {
        std::cout << "[plot_score_systematic_envelopes] skipped '" << output_stem
                  << "' because no universes were accumulated.\n";
        return;
    }

    gStyle->SetOptStat(0);

    auto hcv = make_histogram("hcv_" + canvas_name, s.nbins, s.xmin, s.xmax, s.cv);
    hcv->SetLineColor(line_color);
    hcv->SetLineWidth(3);
    hcv->SetFillStyle(0);

    auto g_rms = make_rms_band(s, line_color);
    auto g_env = make_env_band(s, line_color);
    auto g_mean = make_ratio_curve(s, s.mean_ratio, line_color, 1, 3);
    auto g_min = make_ratio_curve(s, s.min_ratio, line_color, 2, 2);
    auto g_max = make_ratio_curve(s, s.max_ratio, line_color, 2, 2);

    TCanvas c(canvas_name.c_str(), canvas_name.c_str(), 980, 900);
    c.SetTopMargin(0.04);
    c.SetLeftMargin(0.0);
    c.SetRightMargin(0.0);
    c.SetBottomMargin(0.0);

    TPad pad_top((canvas_name + "_top").c_str(), "", 0.0, 0.33, 1.0, 1.0);
    TPad pad_bot((canvas_name + "_bot").c_str(), "", 0.0, 0.00, 1.0, 0.33);

    pad_top.SetLeftMargin(0.12);
    pad_top.SetRightMargin(0.04);
    pad_top.SetBottomMargin(0.02);
    if (logy_top)
        pad_top.SetLogy();

    pad_bot.SetLeftMargin(0.12);
    pad_bot.SetRightMargin(0.04);
    pad_bot.SetTopMargin(0.02);
    pad_bot.SetBottomMargin(0.30);
    pad_bot.SetGridy();

    pad_top.Draw();
    pad_bot.Draw();

    pad_top.cd();

    hcv->GetXaxis()->SetLabelSize(0.0);
    hcv->GetXaxis()->SetTitleSize(0.0);
    hcv->GetYaxis()->SetTitle(y_title_top.c_str());
    hcv->GetYaxis()->SetTitleSize(0.05);
    hcv->GetYaxis()->SetLabelSize(0.045);
    hcv->GetYaxis()->SetTitleOffset(1.15);
    hcv->SetTitle("");

    if (logy_top)
    {
        const double ymax = (s.cv_max > 0.0) ? s.cv_max : 1.0;
        const double ymin = std::isfinite(s.cv_min_positive) ? 0.5 * s.cv_min_positive : 1.0e-6;
        hcv->SetMinimum(std::max(1.0e-9, ymin));
        hcv->SetMaximum(std::max(10.0 * ymax, 10.0 * hcv->GetMinimum()));
    }
    else
    {
        hcv->SetMinimum(0.0);
        hcv->SetMaximum((s.cv_max > 0.0) ? 1.35 * s.cv_max : 1.0);
    }

    hcv->Draw("HIST");

    TLegend leg(0.14, 0.72, 0.72, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    std::ostringstream cv_ss;
    cv_ss << "CV (" << s.nuniv << " universes)";
    leg.AddEntry(hcv.get(), cv_ss.str().c_str(), "l");
    if (envelope_mode == "rms")
        leg.AddEntry(g_rms.get(), "Universe/CV RMS band", "f");
    else if (envelope_mode == "minmax")
        leg.AddEntry(g_env.get(), "Universe/CV min-max envelope", "f");
    else
    {
        leg.AddEntry(g_rms.get(), "Universe/CV RMS band", "f");
        leg.AddEntry(g_min.get(), "Universe/CV min-max envelope", "l");
    }
    leg.Draw();

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.12, 0.93, title_text.c_str());

    pad_bot.cd();

    double rlo = s.ratio_min_global;
    double rhi = s.ratio_max_global;
    const double rspan = std::max(0.08, rhi - rlo);
    rlo -= 0.15 * rspan;
    rhi += 0.15 * rspan;
    rlo = std::max(0.0, rlo);
    if (rlo > 0.95)
        rlo = 0.95;
    if (rhi < 1.05)
        rhi = 1.05;

    TH1D frame(("frame_" + canvas_name).c_str(), "", s.nbins, s.xmin, s.xmax);
    frame.SetDirectory(nullptr);
    frame.SetTitle("");
    frame.SetMinimum(rlo);
    frame.SetMaximum(rhi);
    frame.GetXaxis()->SetTitle(x_title.c_str());
    frame.GetYaxis()->SetTitle("Universe / CV");
    frame.GetXaxis()->SetTitleSize(0.11);
    frame.GetXaxis()->SetLabelSize(0.10);
    frame.GetYaxis()->SetTitleSize(0.10);
    frame.GetYaxis()->SetLabelSize(0.09);
    frame.GetXaxis()->SetTitleOffset(1.05);
    frame.GetYaxis()->SetTitleOffset(0.55);
    frame.GetYaxis()->SetNdivisions(505);
    frame.Draw("AXIS");

    TLine line_one(s.xmin, 1.0, s.xmax, 1.0);
    line_one.SetLineStyle(2);
    line_one.SetLineWidth(2);
    line_one.Draw("SAME");

    if (envelope_mode == "minmax")
    {
        g_env->Draw("E2 SAME");
        g_mean->Draw("L SAME");
    }
    else if (envelope_mode == "rms")
    {
        g_rms->Draw("E2 SAME");
        g_mean->Draw("L SAME");
    }
    else
    {
        g_rms->Draw("E2 SAME");
        g_min->Draw("L SAME");
        g_max->Draw("L SAME");
        g_mean->Draw("L SAME");
    }

    frame.Draw("AXIS SAME");

    c.RedrawAxis();

    const auto out = plot_output_file(output_stem).string();
    c.SaveAs(out.c_str());
    std::cout << "[plot_score_systematic_envelopes] saved: " << out << "\n";
}

std::string family_cv_expression_or_default(ROOT::RDF::RNode node,
                                            const std::string &expr)
{
    if (expr.empty())
        return "1.0";

    if (is_simple_identifier(expr) && !has_column(node, expr))
        return "1.0";

    return expr;
}

void print_family_summary(const FamilyConfig &fam,
                          const EnvelopeAccumulator &acc,
                          bool include_ext_in_total_bkg)
{
    const auto total_bkg_cv = acc.total_bkg_cv(include_ext_in_total_bkg);

    const double sig_total = sum_contents(acc.sig_cv);
    const double bkg_total = sum_contents(total_bkg_cv);
    const double ext_total = sum_contents(acc.bkg_ext_cv);

    double max_sig_rms = 0.0;
    double max_bkg_rms = 0.0;

    const auto sig_summary = summarise_envelope(acc.sig_cv, acc.sig_univ, acc.xmin, acc.xmax);
    const auto bkg_summary = summarise_envelope(total_bkg_cv, acc.total_bkg_univ(include_ext_in_total_bkg), acc.xmin, acc.xmax);

    for (double x : sig_summary.rms_ratio)
        if (std::isfinite(x))
            max_sig_rms = std::max(max_sig_rms, x);

    for (double x : bkg_summary.rms_ratio)
        if (std::isfinite(x))
            max_bkg_rms = std::max(max_bkg_rms, x);

    std::cout << "[plot_score_systematic_envelopes] family=" << fam.key
              << " branch=" << fam.branch
              << " universes=" << acc.nuniv
              << " signal_yield=" << sig_total
              << " total_bkg_yield=" << bkg_total
              << " ext_fixed_yield=" << ext_total
              << " max_signal_bin_rms=" << max_sig_rms
              << " max_bkg_bin_rms=" << max_bkg_rms
              << "\n";
}

} // namespace

int plotFirstInferenceScoreSystematicEnvelopes(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &signal_sel = "is_signal",
    const std::string &base_weight = "w_nominal",
    int nbins = 50,
    double xmin = -15.0,
    double xmax = 15.0,
    bool include_ext_in_total_bkg = true,
    const std::string &envelope_mode = "both",
    const std::string &ppfx_cv_expr = "ppfx_cv",
    const std::string &output_stem = "first_inference_score_syst_envelopes")
{
    return heron::macro::run_with_guard("plotFirstInferenceScoreSystematicEnvelopes", [&]() -> int {
        if (nbins < 1)
            nbins = 1;
        if (xmax < xmin)
            std::swap(xmin, xmax);
        if (xmax == xmin)
            xmax = xmin + 1.0;

        TH1::AddDirectory(false);

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_score_systematic_envelopes] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_score_systematic_envelopes] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode rdf = SelectionService::decorate(el.rdf());

        if (!has_column(rdf, "inf_scores"))
        {
            std::cerr << "[plot_score_systematic_envelopes] missing required branch: inf_scores\n";
            return 1;
        }

        if (is_simple_identifier(base_weight) && !has_column(rdf, base_weight))
        {
            std::cerr << "[plot_score_systematic_envelopes] missing weight column: " << base_weight << "\n";
            return 1;
        }

        if (is_simple_identifier(signal_sel) && !has_column(rdf, signal_sel))
        {
            std::cerr << "[plot_score_systematic_envelopes] missing signal column/expression: " << signal_sel << "\n";
            return 1;
        }

        ROOT::RDF::RNode base = rdf.Define(
                                     "__score__",
                                     [](const ROOT::RVec<float> &scores) {
                                         return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                     },
                                     {"inf_scores"})
                                   .Define("__w__", base_weight);

        if (!base_sel.empty())
        {
            if (is_simple_identifier(base_sel) && !has_column(rdf, base_sel))
            {
                std::cerr << "[plot_score_systematic_envelopes] selection column '" << base_sel
                          << "' is missing; skipping extra selection.\n";
            }
            else
            {
                base = base.Filter(base_sel);
            }
        }

        base = base.Filter(
            [](double s) {
                return std::isfinite(s);
            },
            {"__score__"});

        auto mask_ext = el.mask_for_ext();
        auto mask_mc_like = el.mask_for_mc_like();

        ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext, "sample_id");

        ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc_like, "sample_id")
                                       .Filter(
                                           [mask_ext](int sid) {
                                               return !(sid >= 0 &&
                                                        sid < static_cast<int>(mask_ext->size()) &&
                                                        (*mask_ext)[static_cast<std::size_t>(sid)]);
                                           },
                                           {"sample_id"})
                                       .Define("__is_signal__", signal_sel);

        std::string mode = envelope_mode;
        std::transform(mode.begin(), mode.end(), mode.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (mode != "rms" && mode != "minmax" && mode != "both")
            mode = "both";

        std::vector<FamilyConfig> families{
            {"genie", "GENIE", "weightsGenie", ""},
            {"flux", "Flux", "weightsFlux", ""},
            {"reint", "Reinteraction", "weightsReint", ""},
            {"ppfx", "PPFX", "weightsPPFX", ppfx_cv_expr},
        };

        for (const auto &fam : families)
        {
            if (!has_column(node_mc, fam.branch))
            {
                std::cout << "[plot_score_systematic_envelopes] skipping family '" << fam.key
                          << "' because branch '" << fam.branch << "' is absent.\n";
                continue;
            }

            const std::string cv_expr = family_cv_expression_or_default(node_mc, fam.cv_expr);

            ROOT::RDF::RNode node_mc_fam = node_mc.Define("__family_cv__", cv_expr);

            EnvelopeAccumulator acc(nbins, xmin, xmax);
            std::mutex acc_mutex;

            node_mc_fam.Foreach(
                [&acc, &acc_mutex, &fam](double score,
                                         double w_base,
                                         bool is_signal,
                                         double family_cv,
                                         const ROOT::RVec<unsigned short> &uvec) {
                    std::lock_guard<std::mutex> lock(acc_mutex);
                    acc.fill_mc(is_signal, score, w_base, family_cv, uvec, fam.key);
                },
                {"__score__", "__w__", "__is_signal__", "__family_cv__", fam.branch});

            if (include_ext_in_total_bkg)
            {
                node_ext.Foreach(
                    [&acc, &acc_mutex](double score, double w_base) {
                        std::lock_guard<std::mutex> lock(acc_mutex);
                        acc.fill_ext(score, w_base);
                    },
                    {"__score__", "__w__"});
            }

            if (acc.nuniv == 0)
            {
                std::cout << "[plot_score_systematic_envelopes] skipping family '" << fam.key
                          << "' because no universe entries were found in branch '" << fam.branch << "'.\n";
                continue;
            }

            print_family_summary(fam, acc, include_ext_in_total_bkg);

            const auto total_bkg_cv = acc.total_bkg_cv(include_ext_in_total_bkg);
            const auto total_bkg_univ = acc.total_bkg_univ(include_ext_in_total_bkg);

            const auto sig_tail_cv = reverse_cumulative(acc.sig_cv);
            const auto sig_tail_univ = reverse_cumulative_set(acc.sig_univ);

            const auto bkg_tail_cv = reverse_cumulative(total_bkg_cv);
            const auto bkg_tail_univ = reverse_cumulative_set(total_bkg_univ);

            const auto sig_score_summary = summarise_envelope(acc.sig_cv, acc.sig_univ, xmin, xmax);
            const auto bkg_score_summary = summarise_envelope(total_bkg_cv, total_bkg_univ, xmin, xmax);
            const auto sig_tail_summary = summarise_envelope(sig_tail_cv, sig_tail_univ, xmin, xmax);
            const auto bkg_tail_summary = summarise_envelope(bkg_tail_cv, bkg_tail_univ, xmin, xmax);

            const std::string sig_title = fam.title + " envelope: signal score histogram";
            const std::string bkg_title = include_ext_in_total_bkg
                                              ? fam.title + " envelope: total background score histogram (MC varied, EXT fixed)"
                                              : fam.title + " envelope: total background score histogram";
            const std::string sig_tail_title = fam.title + " envelope: signal cumulative tail N(score > t)";
            const std::string bkg_tail_title = include_ext_in_total_bkg
                                                   ? fam.title + " envelope: total background cumulative tail N(score > t) (MC varied, EXT fixed)"
                                                   : fam.title + " envelope: total background cumulative tail N(score > t)";

            draw_single_envelope_plot(sig_score_summary,
                                      "c_" + fam.key + "_score_signal",
                                      sig_title,
                                      "Inference score [0]",
                                      "Weighted events / bin",
                                      output_stem + "_" + fam.key + "_score_signal",
                                      kBlue + 1,
                                      true,
                                      mode);

            draw_single_envelope_plot(bkg_score_summary,
                                      "c_" + fam.key + "_score_total_bkg",
                                      bkg_title,
                                      "Inference score [0]",
                                      "Weighted events / bin",
                                      output_stem + "_" + fam.key + "_score_total_bkg",
                                      kRed + 1,
                                      true,
                                      mode);

            draw_single_envelope_plot(sig_tail_summary,
                                      "c_" + fam.key + "_tail_signal",
                                      sig_tail_title,
                                      "Inference score threshold t",
                                      "Weighted events with score > t",
                                      output_stem + "_" + fam.key + "_tail_signal",
                                      kBlue + 1,
                                      true,
                                      mode);

            draw_single_envelope_plot(bkg_tail_summary,
                                      "c_" + fam.key + "_tail_total_bkg",
                                      bkg_tail_title,
                                      "Inference score threshold t",
                                      "Weighted events with score > t",
                                      output_stem + "_" + fam.key + "_tail_total_bkg",
                                      kRed + 1,
                                      true,
                                      mode);
        }

        std::cout << "[plot_score_systematic_envelopes] done\n";
        return 0;
    });
}

int plot_score_systematic_envelopes(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &signal_sel = "is_signal",
    const std::string &base_weight = "w_nominal",
    int nbins = 50,
    double xmin = -15.0,
    double xmax = 15.0,
    bool include_ext_in_total_bkg = true,
    const std::string &envelope_mode = "both",
    const std::string &ppfx_cv_expr = "ppfx_cv",
    const std::string &output_stem = "plot_score_systematic_envelopes")
{
    return plotFirstInferenceScoreSystematicEnvelopes(event_list_path,
                                                      base_sel,
                                                      signal_sel,
                                                      base_weight,
                                                      nbins,
                                                      xmin,
                                                      xmax,
                                                      include_ext_in_total_bkg,
                                                      envelope_mode,
                                                      ppfx_cv_expr,
                                                      output_stem);
}
