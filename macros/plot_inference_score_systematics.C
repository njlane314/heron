#ifndef PLOT_INFERENCE_SCORE_SYSTEMATICS_CXX
#define PLOT_INFERENCE_SCORE_SYSTEMATICS_CXX

// macros/plot_inference_score_systematics.C
//
// Systematics suite for the principal evidential axis inf_scores[0].
//
// Produces, for each available multisim source:
//   1) rate+shape fractional uncertainty vs score
//   2) shape-only fractional uncertainty vs score
//   3) cumulative-tail fractional uncertainty vs threshold
//   4) covariance matrix on the score axis
//   5) correlation matrix on the score axis
//   6) shape-only covariance matrix on the score axis
//   7) shape-only correlation matrix on the score axis
//   8) stability of log p(score|S) - log p(score|B) across universes
//
// And, if weightsGenieUp / weightsGenieDn are present:
//   9) one response plot per GENIE knob
//
// Run with:
//   ./heron macro plot_inference_score_systematics.C
//
// Or:
//   ./heron macro plot_inference_score_systematics.C \
//     'plot_inference_score_systematics("./scratch/out/event_list_myana.root")'
//
// Example:
//   ./heron macro plot_inference_score_systematics.C \
//     'plot_inference_score_systematics("./scratch/out/event_list_myana.root", "sel_muon", "is_signal", "w_nominal", 60, -15, 15, "score_systs", true)'

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
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
#include "Plotter.hh"
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

std::string slugify(std::string s)
{
    for (char &c : s)
    {
        if (!std::isalnum(static_cast<unsigned char>(c)))
            c = '_';
        else
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    while (!s.empty() && s.back() == '_')
        s.pop_back();
    while (!s.empty() && s.front() == '_')
        s.erase(s.begin());
    return s.empty() ? "plot" : s;
}

double decode_multisim(unsigned short w)
{
    return 0.001 * static_cast<double>(w);
}

double decode_multisim(double w)
{
    if (!std::isfinite(w))
        return 0.0;
    return std::max(0.0, w);
}

int score_bin(double s, int nbins, double xmin, double xmax)
{
    if (!std::isfinite(s))
        return -1;
    if (s < xmin || s >= xmax)
        return -1;
    const double x = (s - xmin) / (xmax - xmin);
    int b = static_cast<int>(std::floor(x * nbins));
    if (b < 0)
        b = 0;
    if (b >= nbins)
        b = nbins - 1;
    return b;
}

double bin_center(int b, int nbins, double xmin, double xmax)
{
    return xmin + (static_cast<double>(b) + 0.5) * (xmax - xmin) / static_cast<double>(nbins);
}

double vec_sum(const std::vector<double> &v)
{
    return std::accumulate(v.begin(), v.end(), 0.0);
}

std::vector<double> cumulative_tail(const std::vector<double> &v)
{
    std::vector<double> out(v.size(), 0.0);
    double acc = 0.0;
    for (int i = static_cast<int>(v.size()) - 1; i >= 0; --i)
    {
        acc += v[static_cast<std::size_t>(i)];
        out[static_cast<std::size_t>(i)] = acc;
    }
    return out;
}

std::vector<double> extract_universe_hist(const std::vector<double> &flat, int u, int nbins)
{
    std::vector<double> out(static_cast<std::size_t>(nbins), 0.0);
    const std::size_t off = static_cast<std::size_t>(u) * static_cast<std::size_t>(nbins);
    for (int b = 0; b < nbins; ++b)
        out[static_cast<std::size_t>(b)] = flat[off + static_cast<std::size_t>(b)];
    return out;
}

struct SpreadSummary
{
    std::vector<double> cv;
    std::vector<double> mean;
    std::vector<double> rms;
    std::vector<double> min;
    std::vector<double> max;
    std::vector<double> frac;
};

SpreadSummary compute_spread(const std::vector<double> &cv_raw,
                             const std::vector<double> &flat_raw,
                             int nuniv,
                             int nbins,
                             bool shape_only,
                             bool use_tail)
{
    SpreadSummary out;
    out.cv = use_tail ? cumulative_tail(cv_raw) : cv_raw;
    if (shape_only)
    {
        const double norm = vec_sum(out.cv);
        if (norm > 0.0)
            for (double &x : out.cv)
                x /= norm;
        else
            for (double &x : out.cv)
                x = 0.0;
    }

    out.mean.assign(static_cast<std::size_t>(nbins), 0.0);
    out.rms.assign(static_cast<std::size_t>(nbins), 0.0);
    out.min.assign(static_cast<std::size_t>(nbins), std::numeric_limits<double>::infinity());
    out.max.assign(static_cast<std::size_t>(nbins), -std::numeric_limits<double>::infinity());
    out.frac.assign(static_cast<std::size_t>(nbins), 0.0);

    std::vector<double> sum(static_cast<std::size_t>(nbins), 0.0);
    std::vector<double> sum2(static_cast<std::size_t>(nbins), 0.0);
    std::vector<int> count(static_cast<std::size_t>(nbins), 0);

    for (int u = 0; u < nuniv; ++u)
    {
        std::vector<double> vals = extract_universe_hist(flat_raw, u, nbins);
        if (use_tail)
            vals = cumulative_tail(vals);

        if (shape_only)
        {
            const double norm = vec_sum(vals);
            if (!(norm > 0.0))
                continue;
            for (double &x : vals)
                x /= norm;
        }

        for (int b = 0; b < nbins; ++b)
        {
            const double x = vals[static_cast<std::size_t>(b)];
            if (!std::isfinite(x))
                continue;
            sum[static_cast<std::size_t>(b)] += x;
            sum2[static_cast<std::size_t>(b)] += x * x;
            out.min[static_cast<std::size_t>(b)] = std::min(out.min[static_cast<std::size_t>(b)], x);
            out.max[static_cast<std::size_t>(b)] = std::max(out.max[static_cast<std::size_t>(b)], x);
            count[static_cast<std::size_t>(b)] += 1;
        }
    }

    for (int b = 0; b < nbins; ++b)
    {
        const int n = count[static_cast<std::size_t>(b)];
        if (n > 0)
        {
            const double mean = sum[static_cast<std::size_t>(b)] / static_cast<double>(n);
            const double var = std::max(0.0, sum2[static_cast<std::size_t>(b)] / static_cast<double>(n) - mean * mean);
            out.mean[static_cast<std::size_t>(b)] = mean;
            out.rms[static_cast<std::size_t>(b)] = std::sqrt(var);
        }
        else
        {
            out.mean[static_cast<std::size_t>(b)] = 0.0;
            out.rms[static_cast<std::size_t>(b)] = 0.0;
            out.min[static_cast<std::size_t>(b)] = 0.0;
            out.max[static_cast<std::size_t>(b)] = 0.0;
        }

        const double cv = out.cv[static_cast<std::size_t>(b)];
        out.frac[static_cast<std::size_t>(b)] = (cv > 0.0) ? (out.rms[static_cast<std::size_t>(b)] / cv) : 0.0;
    }

    return out;
}

struct MatrixSummary
{
    std::vector<double> cov;
    std::vector<double> corr;
    int nvalid = 0;
};

MatrixSummary compute_covcorr(const std::vector<double> &flat_raw,
                              int nuniv,
                              int nbins,
                              bool shape_only,
                              bool use_tail)
{
    MatrixSummary out;
    out.cov.assign(static_cast<std::size_t>(nbins * nbins), 0.0);
    out.corr.assign(static_cast<std::size_t>(nbins * nbins), 0.0);

    std::vector<std::vector<double>> universes;
    universes.reserve(static_cast<std::size_t>(nuniv));

    for (int u = 0; u < nuniv; ++u)
    {
        std::vector<double> vals = extract_universe_hist(flat_raw, u, nbins);
        if (use_tail)
            vals = cumulative_tail(vals);

        if (shape_only)
        {
            const double norm = vec_sum(vals);
            if (!(norm > 0.0))
                continue;
            for (double &x : vals)
                x /= norm;
        }

        universes.push_back(std::move(vals));
    }

    out.nvalid = static_cast<int>(universes.size());
    if (out.nvalid < 2)
        return out;

    std::vector<double> mean(static_cast<std::size_t>(nbins), 0.0);
    for (const auto &u : universes)
        for (int b = 0; b < nbins; ++b)
            mean[static_cast<std::size_t>(b)] += u[static_cast<std::size_t>(b)];
    for (double &x : mean)
        x /= static_cast<double>(out.nvalid);

    for (const auto &u : universes)
    {
        for (int i = 0; i < nbins; ++i)
        {
            const double di = u[static_cast<std::size_t>(i)] - mean[static_cast<std::size_t>(i)];
            for (int j = 0; j < nbins; ++j)
            {
                const double dj = u[static_cast<std::size_t>(j)] - mean[static_cast<std::size_t>(j)];
                out.cov[static_cast<std::size_t>(i * nbins + j)] += di * dj;
            }
        }
    }

    for (double &x : out.cov)
        x /= static_cast<double>(out.nvalid - 1);

    for (int i = 0; i < nbins; ++i)
    {
        const double vii = out.cov[static_cast<std::size_t>(i * nbins + i)];
        for (int j = 0; j < nbins; ++j)
        {
            const double vjj = out.cov[static_cast<std::size_t>(j * nbins + j)];
            const double denom = std::sqrt(std::max(0.0, vii * vjj));
            out.corr[static_cast<std::size_t>(i * nbins + j)] =
                (denom > 0.0) ? (out.cov[static_cast<std::size_t>(i * nbins + j)] / denom) : 0.0;
        }
    }

    return out;
}

struct LLRSpread
{
    std::vector<double> x;
    std::vector<double> y_cv;
    std::vector<double> y_mean;
    std::vector<double> y_rms;
};

LLRSpread compute_llr_spread(const std::vector<double> &cv_sig,
                             const std::vector<double> &cv_bkg,
                             const std::vector<double> &flat_sig,
                             const std::vector<double> &flat_bkg,
                             int nuniv,
                             int nbins,
                             double xmin,
                             double xmax)
{
    LLRSpread out;
    out.x.resize(static_cast<std::size_t>(nbins), 0.0);
    out.y_cv.resize(static_cast<std::size_t>(nbins), std::numeric_limits<double>::quiet_NaN());
    out.y_mean.resize(static_cast<std::size_t>(nbins), std::numeric_limits<double>::quiet_NaN());
    out.y_rms.resize(static_cast<std::size_t>(nbins), 0.0);

    const double cv_s = vec_sum(cv_sig);
    const double cv_b = vec_sum(cv_bkg);

    for (int b = 0; b < nbins; ++b)
    {
        out.x[static_cast<std::size_t>(b)] = bin_center(b, nbins, xmin, xmax);
        const double ns = cv_sig[static_cast<std::size_t>(b)];
        const double nb = cv_bkg[static_cast<std::size_t>(b)];
        if (cv_s > 0.0 && cv_b > 0.0 && ns > 0.0 && nb > 0.0)
            out.y_cv[static_cast<std::size_t>(b)] = std::log((ns / cv_s) / (nb / cv_b));
    }

    std::vector<double> sum(static_cast<std::size_t>(nbins), 0.0);
    std::vector<double> sum2(static_cast<std::size_t>(nbins), 0.0);
    std::vector<int> count(static_cast<std::size_t>(nbins), 0);

    for (int u = 0; u < nuniv; ++u)
    {
        const std::vector<double> hs = extract_universe_hist(flat_sig, u, nbins);
        const std::vector<double> hb = extract_universe_hist(flat_bkg, u, nbins);

        const double ts = vec_sum(hs);
        const double tb = vec_sum(hb);
        if (!(ts > 0.0) || !(tb > 0.0))
            continue;

        for (int b = 0; b < nbins; ++b)
        {
            const double ns = hs[static_cast<std::size_t>(b)];
            const double nb = hb[static_cast<std::size_t>(b)];
            if (!(ns > 0.0) || !(nb > 0.0))
                continue;
            const double y = std::log((ns / ts) / (nb / tb));
            sum[static_cast<std::size_t>(b)] += y;
            sum2[static_cast<std::size_t>(b)] += y * y;
            count[static_cast<std::size_t>(b)] += 1;
        }
    }

    for (int b = 0; b < nbins; ++b)
    {
        const int n = count[static_cast<std::size_t>(b)];
        if (n > 0)
        {
            const double mean = sum[static_cast<std::size_t>(b)] / static_cast<double>(n);
            const double var = std::max(0.0, sum2[static_cast<std::size_t>(b)] / static_cast<double>(n) - mean * mean);
            out.y_mean[static_cast<std::size_t>(b)] = mean;
            out.y_rms[static_cast<std::size_t>(b)] = std::sqrt(var);
        }
    }

    return out;
}

std::unique_ptr<TH1D> make_hist(const std::string &name,
                                const std::string &title,
                                const std::vector<double> &vals,
                                int nbins,
                                double xmin,
                                double xmax)
{
    auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), nbins, xmin, xmax);
    h->SetDirectory(nullptr);
    for (int b = 0; b < nbins; ++b)
        h->SetBinContent(b + 1, vals[static_cast<std::size_t>(b)]);
    return h;
}

void style_hist(TH1D &h, int color, int style = 1)
{
    h.SetLineColor(color);
    h.SetMarkerColor(color);
    h.SetLineWidth(3);
    h.SetLineStyle(style);
    h.SetMarkerStyle(0);
}

void save_triplet_plot(const std::string &canvas_name,
                       const std::string &title,
                       const std::string &xlabel,
                       const std::string &ylabel,
                       const std::vector<double> &sig,
                       const std::vector<double> &bkg,
                       const std::vector<double> &tot,
                       int nbins,
                       double xmin,
                       double xmax,
                       const std::string &out_stem)
{
    auto hs = make_hist("h_sig_" + canvas_name, ";" + xlabel + ";" + ylabel, sig, nbins, xmin, xmax);
    auto hb = make_hist("h_bkg_" + canvas_name, ";" + xlabel + ";" + ylabel, bkg, nbins, xmin, xmax);
    auto ht = make_hist("h_tot_" + canvas_name, ";" + xlabel + ";" + ylabel, tot, nbins, xmin, xmax);

    style_hist(*hs, kRed + 1, 1);
    style_hist(*hb, kBlue + 1, 2);
    style_hist(*ht, kBlack, 1);

    double ymax = 0.0;
    for (double x : sig)
        ymax = std::max(ymax, x);
    for (double x : bkg)
        ymax = std::max(ymax, x);
    for (double x : tot)
        ymax = std::max(ymax, x);

    TCanvas c(canvas_name.c_str(), title.c_str(), 1100, 800);
    c.SetLeftMargin(0.12);
    c.SetBottomMargin(0.12);

    hs->SetMinimum(0.0);
    hs->SetMaximum((ymax > 0.0) ? 1.20 * ymax : 1.0);
    hs->SetTitle(title.c_str());
    hs->Draw("HIST");
    hb->Draw("HIST SAME");
    ht->Draw("HIST SAME");

    TLegend leg(0.62, 0.73, 0.90, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hs.get(), "signal", "l");
    leg.AddEntry(hb.get(), "background", "l");
    leg.AddEntry(ht.get(), "total MC", "l");
    leg.Draw();

    c.RedrawAxis();
    const auto out = plot_output_file(out_stem).string();
    c.SaveAs(out.c_str());
    std::cout << "[plot_inference_score_systematics] saved: " << out << "\n";
}

void save_matrix_plot(const std::string &canvas_name,
                      const std::string &title,
                      const std::vector<double> &mat,
                      int nbins,
                      double xmin,
                      double xmax,
                      const std::string &zlabel,
                      const std::string &out_stem,
                      bool force_corr_range)
{
    TCanvas c(canvas_name.c_str(), title.c_str(), 1000, 850);
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.16);
    c.SetBottomMargin(0.12);

    TH2D h(("h2_" + canvas_name).c_str(),
           (";" + std::string("Inference score [0] bin centre") + ";" + std::string("Inference score [0] bin centre")).c_str(),
           nbins,
           xmin,
           xmax,
           nbins,
           xmin,
           xmax);

    for (int i = 0; i < nbins; ++i)
        for (int j = 0; j < nbins; ++j)
            h.SetBinContent(i + 1, j + 1, mat[static_cast<std::size_t>(i * nbins + j)]);

    h.SetTitle(title.c_str());
    if (force_corr_range)
    {
        h.SetMinimum(-1.0);
        h.SetMaximum(+1.0);
    }
    h.GetZaxis()->SetTitle(zlabel.c_str());
    h.Draw("COLZ");

    c.RedrawAxis();
    const auto out = plot_output_file(out_stem).string();
    c.SaveAs(out.c_str());
    std::cout << "[plot_inference_score_systematics] saved: " << out << "\n";
}

void save_llr_stability_plot(const std::string &canvas_name,
                             const std::string &title,
                             const LLRSpread &llr,
                             double xmin,
                             double xmax,
                             const std::string &out_stem)
{
    std::vector<double> xb, yb, exb, eyb;
    std::vector<double> xc, yc;
    std::vector<double> xm, ym;

    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < llr.x.size(); ++i)
    {
        const double x = llr.x[i];
        const double ycv = llr.y_cv[i];
        const double ymn = llr.y_mean[i];
        const double yr = llr.y_rms[i];

        if (std::isfinite(ymn))
        {
            xb.push_back(x);
            yb.push_back(ymn);
            exb.push_back(0.0);
            eyb.push_back(yr);
            xm.push_back(x);
            ym.push_back(ymn);
            ymin = std::min(ymin, ymn - yr);
            ymax = std::max(ymax, ymn + yr);
        }

        if (std::isfinite(ycv))
        {
            xc.push_back(x);
            yc.push_back(ycv);
            ymin = std::min(ymin, ycv);
            ymax = std::max(ymax, ycv);
        }
    }

    if (xb.empty() && xc.empty())
        return;

    if (!std::isfinite(ymin) || !std::isfinite(ymax))
    {
        ymin = -5.0;
        ymax = 5.0;
    }

    const double pad = 0.12 * std::max(1.0, ymax - ymin);

    TCanvas c(canvas_name.c_str(), title.c_str(), 1100, 800);
    c.SetLeftMargin(0.12);
    c.SetBottomMargin(0.12);

    TH1D frame(("frame_" + canvas_name).c_str(),
               ";Inference score [0];log p(score|S) - log p(score|B)",
               100,
               xmin,
               xmax);
    frame.SetTitle(title.c_str());
    frame.SetMinimum(ymin - pad);
    frame.SetMaximum(ymax + pad);
    frame.Draw("AXIS");

    if (!xb.empty())
    {
        TGraphErrors gband(static_cast<int>(xb.size()), xb.data(), yb.data(), exb.data(), eyb.data());
        gband.SetFillColorAlpha(kAzure - 9, 0.45);
        gband.SetLineColor(kAzure + 2);
        gband.SetLineWidth(2);
        gband.Draw("3 SAME");

        TGraph gmean(static_cast<int>(xm.size()), xm.data(), ym.data());
        gmean.SetLineColor(kAzure + 3);
        gmean.SetLineWidth(3);
        gmean.Draw("L SAME");
    }

    if (!xc.empty())
    {
        TGraph gcv(static_cast<int>(xc.size()), xc.data(), yc.data());
        gcv.SetLineColor(kBlack);
        gcv.SetLineWidth(3);
        gcv.Draw("L SAME");
    }

    TLegend leg(0.13, 0.74, 0.57, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry((TObject *)0, "CV curve: black", "");
    leg.AddEntry((TObject *)0, "Universe mean #pm RMS: blue band", "");
    leg.Draw();

    c.RedrawAxis();
    const auto out = plot_output_file(out_stem).string();
    c.SaveAs(out.c_str());
    std::cout << "[plot_inference_score_systematics] saved: " << out << "\n";
}

struct EnsembleAccumulator
{
    int nbins = 0;
    int nuniv = 0;
    double xmin = 0.0;
    double xmax = 0.0;

    std::vector<double> cv_sig;
    std::vector<double> cv_bkg;
    std::vector<double> cv_tot;

    std::vector<double> flat_sig;
    std::vector<double> flat_bkg;
    std::vector<double> flat_tot;

    EnsembleAccumulator(int nbin, double xlo, double xhi)
        : nbins(nbin), xmin(xlo), xmax(xhi),
          cv_sig(static_cast<std::size_t>(nbin), 0.0),
          cv_bkg(static_cast<std::size_t>(nbin), 0.0),
          cv_tot(static_cast<std::size_t>(nbin), 0.0)
    {
    }

    void init(int nu)
    {
        nuniv = std::max(0, nu);
        flat_sig.assign(static_cast<std::size_t>(nuniv * nbins), 0.0);
        flat_bkg.assign(static_cast<std::size_t>(nuniv * nbins), 0.0);
        flat_tot.assign(static_cast<std::size_t>(nuniv * nbins), 0.0);
    }

    template <class Vec>
    void fill(double score, bool is_sig, double wcv, const Vec &weights)
    {
        const int b = score_bin(score, nbins, xmin, xmax);
        if (b < 0)
            return;
        if (!std::isfinite(wcv) || wcv < 0.0)
            return;

        cv_tot[static_cast<std::size_t>(b)] += wcv;
        if (is_sig)
            cv_sig[static_cast<std::size_t>(b)] += wcv;
        else
            cv_bkg[static_cast<std::size_t>(b)] += wcv;

        if (nuniv == 0 && !weights.empty())
            init(static_cast<int>(weights.size()));
        if (nuniv == 0)
            return;

        const std::size_t nu = std::min(static_cast<std::size_t>(nuniv), static_cast<std::size_t>(weights.size()));
        for (std::size_t u = 0; u < nu; ++u)
        {
            const double r = decode_multisim(weights[u]);
            const double w = wcv * r;
            const std::size_t idx = u * static_cast<std::size_t>(nbins) + static_cast<std::size_t>(b);
            flat_tot[idx] += w;
            if (is_sig)
                flat_sig[idx] += w;
            else
                flat_bkg[idx] += w;
        }
    }
};

struct KnobAccumulator
{
    int nbins = 0;
    int nknob = 0;
    double xmin = 0.0;
    double xmax = 0.0;

    std::vector<double> cv_tot;
    std::vector<double> flat_up;
    std::vector<double> flat_dn;

    KnobAccumulator(int nbin, double xlo, double xhi)
        : nbins(nbin), xmin(xlo), xmax(xhi),
          cv_tot(static_cast<std::size_t>(nbin), 0.0)
    {
    }

    void init(int nk)
    {
        nknob = std::max(0, nk);
        flat_up.assign(static_cast<std::size_t>(nknob * nbins), 0.0);
        flat_dn.assign(static_cast<std::size_t>(nknob * nbins), 0.0);
    }

    template <class Vec>
    void fill(double score, double wcv, const Vec &up, const Vec &dn)
    {
        const int b = score_bin(score, nbins, xmin, xmax);
        if (b < 0)
            return;
        if (!std::isfinite(wcv) || wcv < 0.0)
            return;

        cv_tot[static_cast<std::size_t>(b)] += wcv;

        const int nk = static_cast<int>(std::min(static_cast<std::size_t>(up.size()), static_cast<std::size_t>(dn.size())));
        if (nknob == 0 && nk > 0)
            init(nk);
        if (nknob == 0)
            return;

        const int nuse = std::min(nknob, nk);
        for (int k = 0; k < nuse; ++k)
        {
            const double wu = wcv * decode_multisim(up[static_cast<std::size_t>(k)]);
            const double wd = wcv * decode_multisim(dn[static_cast<std::size_t>(k)]);
            const std::size_t idx = static_cast<std::size_t>(k * nbins + b);
            flat_up[idx] += wu;
            flat_dn[idx] += wd;
        }
    }
};

std::vector<std::string> genie_knob_names()
{
    return {
        "AGKYpT1pi_UBGenie",
        "AGKYxF1pi_UBGenie",
        "AhtBY_UBGenie",
        "AxFFCCQEshape_UBGenie",
        "BhtBY_UBGenie",
        "CV1uBY_UBGenie",
        "CV2uBY_UBGenie",
        "DecayAngMEC_UBGenie",
        "EtaNCEL_UBGenie",
        "FrAbs_N_UBGenie",
        "FrAbs_pi_UBGenie",
        "FrCEx_N_UBGenie",
        "FrCEx_pi_UBGenie",
        "FrInel_N_UBGenie",
        "FrInel_pi_UBGenie",
        "FrPiProd_N_UBGenie",
        "FrPiProd_pi_UBGenie",
        "FracDelta_CCMEC_UBGenie",
        "FracPN_CCMEC_UBGenie",
        "MFP_N_UBGenie",
        "MFP_pi_UBGenie",
        "MaCCQE_UBGenie",
        "MaCCRES_UBGenie",
        "MaNCEL_UBGenie",
        "MaNCRES_UBGenie",
        "MvCCRES_UBGenie",
        "MvNCRES_UBGenie",
        "NonRESBGvbarnCC1pi_UBGenie",
        "NonRESBGvbarnCC2pi_UBGenie",
        "NonRESBGvbarnNC1pi_UBGenie",
        "NonRESBGvbarnNC2pi_UBGenie",
        "NonRESBGvbarpCC1pi_UBGenie",
        "NonRESBGvbarpCC2pi_UBGenie",
        "NonRESBGvbarpNC1pi_UBGenie",
        "NonRESBGvbarpNC2pi_UBGenie",
        "NonRESBGvnCC1pi_UBGenie",
        "NonRESBGvnCC2pi_UBGenie",
        "NonRESBGvnNC1pi_UBGenie",
        "NonRESBGvnNC2pi_UBGenie",
        "NonRESBGvpCC1pi_UBGenie",
        "NonRESBGvpCC2pi_UBGenie",
        "NonRESBGvpNC1pi_UBGenie",
        "NonRESBGvpNC2pi_UBGenie",
        "NormCCMEC_UBGenie",
        "NormNCMEC_UBGenie",
        "RDecBR1eta_UBGenie",
        "RDecBR1gamma_UBGenie",
        "RPA_CCQE_UBGenie",
        "Theta_Delta2Npi_UBGenie",
        "TunedCentralValue_UBGenie",
        "VecFFCCQEshape_UBGenie",
        "XSecShape_CCMEC_UBGenie",
        "splines_general_Spline"};
}

void save_knob_response_plot(const std::string &name,
                             const std::vector<double> &cv,
                             const std::vector<double> &up,
                             const std::vector<double> &dn,
                             int nbins,
                             double xmin,
                             double xmax,
                             const std::string &out_stem)
{
    std::vector<double> rup(static_cast<std::size_t>(nbins), 0.0);
    std::vector<double> rdn(static_cast<std::size_t>(nbins), 0.0);
    std::vector<double> delt(static_cast<std::size_t>(nbins), 0.0);

    double top_min = std::numeric_limits<double>::infinity();
    double top_max = -std::numeric_limits<double>::infinity();
    double bot_abs = 0.0;

    for (int b = 0; b < nbins; ++b)
    {
        const double c = cv[static_cast<std::size_t>(b)];
        if (c > 0.0)
        {
            rup[static_cast<std::size_t>(b)] = up[static_cast<std::size_t>(b)] / c;
            rdn[static_cast<std::size_t>(b)] = dn[static_cast<std::size_t>(b)] / c;
            delt[static_cast<std::size_t>(b)] = 0.5 * (up[static_cast<std::size_t>(b)] - dn[static_cast<std::size_t>(b)]) / c;
            top_min = std::min(top_min, std::min(rup[static_cast<std::size_t>(b)], rdn[static_cast<std::size_t>(b)]));
            top_max = std::max(top_max, std::max(rup[static_cast<std::size_t>(b)], rdn[static_cast<std::size_t>(b)]));
            bot_abs = std::max(bot_abs, std::fabs(delt[static_cast<std::size_t>(b)]));
        }
    }

    if (!std::isfinite(top_min) || !std::isfinite(top_max))
    {
        top_min = 0.8;
        top_max = 1.2;
    }

    top_min = std::min(top_min, 1.0);
    top_max = std::max(top_max, 1.0);
    const double top_pad = 0.15 * std::max(0.2, top_max - top_min);
    const double bot_lim = std::max(0.02, 1.20 * bot_abs);

    auto hup = make_hist("hup_" + slugify(name), ";Inference score [0];Var/CV", rup, nbins, xmin, xmax);
    auto hdn = make_hist("hdn_" + slugify(name), ";Inference score [0];Var/CV", rdn, nbins, xmin, xmax);
    auto hde = make_hist("hde_" + slugify(name), ";Inference score [0];(Up-Dn)/(2CV)", delt, nbins, xmin, xmax);

    style_hist(*hup, kRed + 1, 1);
    style_hist(*hdn, kBlue + 1, 2);
    style_hist(*hde, kBlack, 1);

    TCanvas c(("c_" + slugify(name)).c_str(), name.c_str(), 1100, 900);

    TPad p1("p1", "", 0.0, 0.33, 1.0, 1.0);
    TPad p2("p2", "", 0.0, 0.00, 1.0, 0.33);

    p1.SetLeftMargin(0.12);
    p1.SetRightMargin(0.05);
    p1.SetBottomMargin(0.02);

    p2.SetLeftMargin(0.12);
    p2.SetRightMargin(0.05);
    p2.SetTopMargin(0.02);
    p2.SetBottomMargin(0.28);

    p1.Draw();
    p2.Draw();

    p1.cd();
    hup->SetTitle(name.c_str());
    hup->SetMinimum(top_min - top_pad);
    hup->SetMaximum(top_max + top_pad);
    hup->Draw("HIST");
    hdn->Draw("HIST SAME");

    TLine l1(xmin, 1.0, xmax, 1.0);
    l1.SetLineStyle(2);
    l1.Draw("SAME");

    TLegend leg(0.64, 0.76, 0.91, 0.90);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hup.get(), "up/CV", "l");
    leg.AddEntry(hdn.get(), "dn/CV", "l");
    leg.Draw();

    p2.cd();
    hde->SetMinimum(-bot_lim);
    hde->SetMaximum(+bot_lim);
    hde->Draw("HIST");

    TLine l0(xmin, 0.0, xmax, 0.0);
    l0.SetLineStyle(2);
    l0.Draw("SAME");

    c.RedrawAxis();
    const auto out = plot_output_file(out_stem).string();
    c.SaveAs(out.c_str());
    std::cout << "[plot_inference_score_systematics] saved: " << out << "\n";
}

} // namespace

int plotInferenceScoreSystematics(const std::string &event_list_path = "",
                                  const std::string &base_sel = "sel_muon",
                                  const std::string &signal_sel = "is_signal",
                                  const std::string &mc_weight = "w_nominal",
                                  int nbins = 60,
                                  double xmin = -15.0,
                                  double xmax = 15.0,
                                  const std::string &output_stem = "inference_score_systematics",
                                  bool make_knob_plots = true)
{
    return heron::macro::run_with_guard("plotInferenceScoreSystematics", [&]() -> int {
        if (nbins < 1)
            nbins = 1;
        if (xmax < xmin)
            std::swap(xmin, xmax);
        if (xmax == xmin)
            xmax = xmin + 1.0;

        TH1::SetDefaultSumw2();
        gStyle->SetOptStat(0);

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_inference_score_systematics] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_inference_score_systematics] input is not an event-list root file: " << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode rdf = SelectionService::decorate(el.rdf())
                                   .Define("__score__",
                                           [](const ROOT::RVec<float> &scores) {
                                               return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                           },
                                           {"inf_scores"});

        if (!has_column(rdf, mc_weight))
        {
            std::cerr << "[plot_inference_score_systematics] missing weight column/expression: '" << mc_weight << "'\n";
            return 1;
        }

        ROOT::RDF::RNode node = filter_by_sample_mask(rdf, el.mask_for_mc_like(), "sample_id")
                                    .Define("__wcv__", mc_weight)
                                    .Define("__is_signal__", signal_sel);

        if (!base_sel.empty())
        {
            if (has_column(node, base_sel))
                node = node.Filter([](bool pass) { return pass; }, {base_sel});
            else
                node = node.Filter(base_sel);
        }

        node = node.Filter(
                   [xmin, xmax](double s) {
                       return std::isfinite(s) && s >= xmin && s < xmax;
                   },
                   {"__score__"})
                   .Filter([](double w) { return std::isfinite(w) && w >= 0.0; }, {"__wcv__"});

        struct SourceSpec
        {
            std::string label;
            std::string branch;
        };

        std::vector<SourceSpec> sources;
        if (has_column(node, "weightsGenie"))
            sources.push_back({"GENIE multisim", "weightsGenie"});
        if (has_column(node, "weightsReint"))
            sources.push_back({"Reinteraction multisim", "weightsReint"});
        if (has_column(node, "weightsPPFX"))
            sources.push_back({"PPFX multisim", "weightsPPFX"});
        if (has_column(node, "weightsFlux"))
            sources.push_back({"Flux multisim", "weightsFlux"});

        if (sources.empty())
        {
            std::cerr << "[plot_inference_score_systematics] no direct multisim vector branches were found.\n";
            std::cerr << "  Expected one or more of: weightsGenie, weightsReint, weightsPPFX, weightsFlux\n";
            std::cerr << "  If your flux universes only live inside the map branch, this macro needs the map fallback added.\n";
            return 1;
        }

        Plotter plotter;
        plotter.set_global_style();

        for (const auto &src : sources)
        {
            const std::string src_slug = slugify(src.label);
            EnsembleAccumulator acc(nbins, xmin, xmax);

            if (src.branch == "weightsGenie")
            {
                node.Foreach(
                    [&acc](double s, bool is_sig, double wcv, const ROOT::RVec<unsigned short> &w) {
                        acc.fill(s, is_sig, wcv, w);
                    },
                    {"__score__", "__is_signal__", "__wcv__", src.branch});
            }
            else if (src.branch == "weightsReint")
            {
                node.Foreach(
                    [&acc](double s, bool is_sig, double wcv, const ROOT::RVec<unsigned short> &w) {
                        acc.fill(s, is_sig, wcv, w);
                    },
                    {"__score__", "__is_signal__", "__wcv__", src.branch});
            }
            else if (src.branch == "weightsPPFX")
            {
                node.Foreach(
                    [&acc](double s, bool is_sig, double wcv, const ROOT::RVec<unsigned short> &w) {
                        acc.fill(s, is_sig, wcv, w);
                    },
                    {"__score__", "__is_signal__", "__wcv__", src.branch});
            }
            else if (src.branch == "weightsFlux")
            {
                node.Foreach(
                    [&acc](double s, bool is_sig, double wcv, const ROOT::RVec<unsigned short> &w) {
                        acc.fill(s, is_sig, wcv, w);
                    },
                    {"__score__", "__is_signal__", "__wcv__", src.branch});
            }

            if (acc.nuniv <= 0)
            {
                std::cout << "[plot_inference_score_systematics] source '" << src.label
                          << "' had no populated universes; skipping.\n";
                continue;
            }

            std::cout << "[plot_inference_score_systematics] source=" << src.label
                      << " universes=" << acc.nuniv << "\n";

            const auto rs_sig = compute_spread(acc.cv_sig, acc.flat_sig, acc.nuniv, nbins, false, false);
            const auto rs_bkg = compute_spread(acc.cv_bkg, acc.flat_bkg, acc.nuniv, nbins, false, false);
            const auto rs_tot = compute_spread(acc.cv_tot, acc.flat_tot, acc.nuniv, nbins, false, false);

            const auto sh_sig = compute_spread(acc.cv_sig, acc.flat_sig, acc.nuniv, nbins, true, false);
            const auto sh_bkg = compute_spread(acc.cv_bkg, acc.flat_bkg, acc.nuniv, nbins, true, false);
            const auto sh_tot = compute_spread(acc.cv_tot, acc.flat_tot, acc.nuniv, nbins, true, false);

            const auto tl_sig = compute_spread(acc.cv_sig, acc.flat_sig, acc.nuniv, nbins, false, true);
            const auto tl_bkg = compute_spread(acc.cv_bkg, acc.flat_bkg, acc.nuniv, nbins, false, true);
            const auto tl_tot = compute_spread(acc.cv_tot, acc.flat_tot, acc.nuniv, nbins, false, true);

            save_triplet_plot("c_" + src_slug + "_rate_shape_frac",
                              src.label + ": rate+shape fractional uncertainty",
                              "Inference score [0]",
                              "RMS / CV",
                              rs_sig.frac,
                              rs_bkg.frac,
                              rs_tot.frac,
                              nbins,
                              xmin,
                              xmax,
                              output_stem + "_" + src_slug + "_rate_shape_fractional_uncertainty");

            save_triplet_plot("c_" + src_slug + "_shape_only_frac",
                              src.label + ": shape-only fractional uncertainty",
                              "Inference score [0]",
                              "RMS(shape) / CV(shape)",
                              sh_sig.frac,
                              sh_bkg.frac,
                              sh_tot.frac,
                              nbins,
                              xmin,
                              xmax,
                              output_stem + "_" + src_slug + "_shape_only_fractional_uncertainty");

            save_triplet_plot("c_" + src_slug + "_tail_frac",
                              src.label + ": cumulative-tail fractional uncertainty",
                              "Threshold t on inference score [0] for N(score > t)",
                              "RMS / CV",
                              tl_sig.frac,
                              tl_bkg.frac,
                              tl_tot.frac,
                              nbins,
                              xmin,
                              xmax,
                              output_stem + "_" + src_slug + "_cumulative_tail_fractional_uncertainty");

            const auto raw_m = compute_covcorr(acc.flat_tot, acc.nuniv, nbins, false, false);
            save_matrix_plot("c_" + src_slug + "_cov_total",
                             src.label + ": covariance matrix on score axis (total MC)",
                             raw_m.cov,
                             nbins,
                             xmin,
                             xmax,
                             "Covariance",
                             output_stem + "_" + src_slug + "_covariance_total",
                             false);

            save_matrix_plot("c_" + src_slug + "_corr_total",
                             src.label + ": correlation matrix on score axis (total MC)",
                             raw_m.corr,
                             nbins,
                             xmin,
                             xmax,
                             "Correlation",
                             output_stem + "_" + src_slug + "_correlation_total",
                             true);

            const auto shp_m = compute_covcorr(acc.flat_tot, acc.nuniv, nbins, true, false);
            save_matrix_plot("c_" + src_slug + "_cov_shape_total",
                             src.label + ": shape-only covariance matrix on score axis (total MC)",
                             shp_m.cov,
                             nbins,
                             xmin,
                             xmax,
                             "Covariance",
                             output_stem + "_" + src_slug + "_shape_only_covariance_total",
                             false);

            save_matrix_plot("c_" + src_slug + "_corr_shape_total",
                             src.label + ": shape-only correlation matrix on score axis (total MC)",
                             shp_m.corr,
                             nbins,
                             xmin,
                             xmax,
                             "Correlation",
                             output_stem + "_" + src_slug + "_shape_only_correlation_total",
                             true);

            const auto llr = compute_llr_spread(acc.cv_sig, acc.cv_bkg, acc.flat_sig, acc.flat_bkg, acc.nuniv, nbins, xmin, xmax);
            save_llr_stability_plot("c_" + src_slug + "_llr_stability",
                                    src.label + ": stability of log p(score|S)-log p(score|B)",
                                    llr,
                                    xmin,
                                    xmax,
                                    output_stem + "_" + src_slug + "_llr_stability");
        }

        if (make_knob_plots && has_column(node, "weightsGenieUp") && has_column(node, "weightsGenieDn"))
        {
            KnobAccumulator kacc(nbins, xmin, xmax);

            node.Foreach(
                [&kacc](double s, double wcv, const ROOT::RVec<unsigned short> &up, const ROOT::RVec<unsigned short> &dn) {
                    kacc.fill(s, wcv, up, dn);
                },
                {"__score__", "__wcv__", "weightsGenieUp", "weightsGenieDn"});

            if (kacc.nknob > 0)
            {
                const auto names = genie_knob_names();
                const int nk = std::min(kacc.nknob, static_cast<int>(names.size()));

                for (int k = 0; k < nk; ++k)
                {
                    std::vector<double> hup(static_cast<std::size_t>(nbins), 0.0);
                    std::vector<double> hdn(static_cast<std::size_t>(nbins), 0.0);

                    for (int b = 0; b < nbins; ++b)
                    {
                        hup[static_cast<std::size_t>(b)] = kacc.flat_up[static_cast<std::size_t>(k * nbins + b)];
                        hdn[static_cast<std::size_t>(b)] = kacc.flat_dn[static_cast<std::size_t>(k * nbins + b)];
                    }

                    save_knob_response_plot(names[static_cast<std::size_t>(k)],
                                            kacc.cv_tot,
                                            hup,
                                            hdn,
                                            nbins,
                                            xmin,
                                            xmax,
                                            output_stem + "_genie_knob_" + slugify(names[static_cast<std::size_t>(k)]));
                }
            }
            else
            {
                std::cout << "[plot_inference_score_systematics] GENIE up/down knob vectors were present but empty.\n";
            }
        }
        else
        {
            std::cout << "[plot_inference_score_systematics] skipping explicit GENIE knob plots.\n";
        }

        std::cout << "[plot_inference_score_systematics] done\n";
        return 0;
    });
}

int plot_inference_score_systematics(const std::string &event_list_path = "",
                                     const std::string &base_sel = "sel_muon",
                                     const std::string &signal_sel = "is_signal",
                                     const std::string &mc_weight = "w_nominal",
                                     int nbins = 60,
                                     double xmin = -15.0,
                                     double xmax = 15.0,
                                     const std::string &output_stem = "inference_score_systematics",
                                     bool make_knob_plots = true)
{
    return plotInferenceScoreSystematics(event_list_path,
                                         base_sel,
                                         signal_sel,
                                         mc_weight,
                                         nbins,
                                         xmin,
                                         xmax,
                                         output_stem,
                                         make_knob_plots);
}

#endif
