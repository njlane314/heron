/* -- C++ -- */
/**
 *  @file  plot/src/AdaptiveBinningService.cpp
 *
 *  @brief Adaptive binning service implementation.
 */

#include "AdaptiveBinningService.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

#include <TAxis.h>
#include <TH1D.h>

#include "PlotDescriptors.hh"

namespace nu
{
namespace
{
constexpr double kEdgeEps = 1e-12;

// ROOT histograms are, by default, registered in the current directory (gDirectory).
// In batch workflows where we own histograms via std::unique_ptr, that implicit ownership
// can lead to double-deletes during teardown (often observed right after TCanvas::Print).
// Guard against that by disabling directory auto-registration in these helpers.
struct AddDirectoryGuard
{
    const Bool_t prev = TH1::AddDirectoryStatus();
    AddDirectoryGuard() { TH1::AddDirectory(kFALSE); }
    ~AddDirectoryGuard() { TH1::AddDirectory(prev); }
};

inline void ensure_sumw2(TH1D &h)
{
    // Avoid ROOT warnings: "Sum of squares of weights structure already created".
    // Also guarantees GetSumw2() is meaningful when using weighted fills.
    if (h.GetSumw2N() == 0)
    {
        h.Sumw2(kTRUE);
    }
}

inline double bin_sumw2(const TH1D &h, int bin)
{
    // Prefer the real stored sumw2 (true for weighted histograms).
    if (h.GetSumw2N() > 0)
    {
        const auto *const p_arr = h.GetSumw2();
        if (p_arr != NULL && bin >= 0 && bin < p_arr->GetSize())
        {
            return p_arr->At(bin);
        }
    }

    // Fallback: assume Poisson (unweighted) => sumw2 ~= sumw.
    // Use abs to be safe if content can be negative.
    return std::abs(h.GetBinContent(bin));
}

inline std::pair<int, int> interior_bins_for_overflow(const TH1D &h, int edge_bins)
{
    const int nb = h.GetNbinsX();
    if (nb <= 0)
    {
        return {1, 1};
    }
    const int n = std::max(0, edge_bins);
    int first = 1 + n;
    int last = nb - n;
    first = std::max(1, std::min(first, nb));
    last = std::max(1, std::min(last, nb));
    if (first > last)
    {
        // Degenerate: if edge_bins is too large, just fold to the edges.
        first = 1;
        last = nb;
    }
    return {first, last};
}

inline void fold_overflow_into(TH1D &h, int first_bin, int last_bin)
{
    const int nb = h.GetNbinsX();
    if (nb <= 0)
    {
        return;
    }
    first_bin = std::max(1, std::min(first_bin, nb));
    last_bin = std::max(1, std::min(last_bin, nb));

    // Underflow -> first_bin
    {
        const double c0 = h.GetBinContent(0);
        const double e0 = h.GetBinError(0);
        const double c1 = h.GetBinContent(first_bin);
        const double e1 = h.GetBinError(first_bin);
        h.SetBinContent(first_bin, c1 + c0);
        h.SetBinError(first_bin, std::hypot(e1, e0));
        h.SetBinContent(0, 0.0);
        h.SetBinError(0, 0.0);
    }

    // Overflow -> last_bin
    {
        const double co = h.GetBinContent(nb + 1);
        const double eo = h.GetBinError(nb + 1);
        const double cn = h.GetBinContent(last_bin);
        const double en = h.GetBinError(last_bin);
        h.SetBinContent(last_bin, cn + co);
        h.SetBinError(last_bin, std::hypot(en, eo));
        h.SetBinContent(nb + 1, 0.0);
        h.SetBinError(nb + 1, 0.0);
    }
}

inline double denom_sumw(double sumw, bool use_abs)
{
    return use_abs ? std::abs(sumw) : sumw;
}

inline bool pass_bin(double sumw, double sumw2, const AdaptiveBinningService::MinStatConfig &cfg)
{
    const bool use_wmin = (cfg.min_sumw > 0.0);
    const bool use_rel = (cfg.max_rel_err > 0.0);

    if (!use_wmin && !use_rel)
    {
        return true;
    }

    const double d = denom_sumw(sumw, cfg.use_abs_sumw);

    if (use_wmin && !(d >= cfg.min_sumw))
    {
        return false;
    }

    if (use_rel)
    {
        if (!(d > 0.0))
        {
            return false;
        }
        const double rel = std::sqrt(std::max(0.0, sumw2)) / d;
        if (!(rel <= cfg.max_rel_err))
        {
            return false;
        }
    }

    return true;
}

inline void log_adaptive_bin_sizes(std::string_view hist_name,
                                   const std::vector<double> &edges)
{
    if (edges.size() < 2)
        return;

    const std::size_t nbins = edges.size() - 1;
    const auto width = [&](std::size_t i) { return edges[i + 1] - edges[i]; };

    double wmin = width(0);
    double wmax = wmin;
    double wsum = wmin;
    for (std::size_t i = 1; i < nbins; ++i)
    {
        const double w = width(i);
        wmin = std::min(wmin, w);
        wmax = std::max(wmax, w);
        wsum += w;
    }
    const double wmean = wsum / static_cast<double>(nbins);

    std::ostringstream msg;
    msg << "[AdaptiveBinningService] Adaptive bins settled for '" << hist_name << "': "
        << nbins << " bins; width min/max/mean = " << wmin << "/" << wmax << "/" << wmean;

    // Avoid dumping thousands of widths to the terminal (hard to read, and can trigger
    // ROOT/Cling instability on some batch nodes).
    constexpr std::size_t kShow = 12;
    msg << "; widths [";
    if (nbins <= 2 * kShow)
    {
        for (std::size_t i = 0; i < nbins; ++i)
        {
            if (i > 0)
                msg << ", ";
            msg << width(i);
        }
    }
    else
    {
        for (std::size_t i = 0; i < kShow; ++i)
        {
            if (i > 0)
                msg << ", ";
            msg << width(i);
        }
        msg << ", ..., ";
        for (std::size_t i = nbins - kShow; i < nbins; ++i)
        {
            if (i > nbins - kShow)
                msg << ", ";
            msg << width(i);
        }
    }
    msg << "]";
    std::clog << msg.str() << '\n';
}

inline std::vector<double> uniform_edges_from(const TH1D &h)
{
    std::vector<double> edges;
    const int nb = h.GetNbinsX();
    if (nb <= 0)
    {
        return edges;
    }

    edges.reserve(static_cast<std::size_t>(nb) + 1);
    const auto *ax = h.GetXaxis();
    edges.push_back(ax->GetBinLowEdge(1));
    for (int i = 1; i <= nb; ++i)
    {
        edges.push_back(ax->GetBinUpEdge(i));
    }
    return edges;
}

} // namespace

AdaptiveBinningService &AdaptiveBinningService::instance()
{
    static AdaptiveBinningService svc;
    return svc;
}

AdaptiveBinningService::MinStatConfig AdaptiveBinningService::config_from(const Options &opt)
{
    MinStatConfig cfg;
    cfg.enabled = opt.adaptive_binning;
    cfg.min_sumw = opt.adaptive_min_sumw;
    cfg.max_rel_err = opt.adaptive_max_relerr;
    cfg.fold_overflow = opt.adaptive_fold_overflow;
    cfg.edge_bins = std::max(0, opt.adaptive_edge_bins);
    cfg.use_abs_sumw = true;
    return cfg;
}

void AdaptiveBinningService::fold_overflow(TH1D &h) const
{
    const int nb = h.GetNbinsX();
    if (nb <= 0)
    {
        return;
    }
    ensure_sumw2(h);
    fold_overflow_into(h, 1, nb);
}

std::unique_ptr<TH1D> AdaptiveBinningService::sum_hists(const std::vector<const TH1D *> &parts,
                                                        std::string_view new_name,
                                                        const MinStatConfig &cfg) const
{
    const TH1D *first = nullptr;
    for (auto *p : parts)
    {
        if (p)
        {
            first = p;
            break;
        }
    }
    if (!first)
    {
        return {};
    }

    const std::string name(new_name);
    AddDirectoryGuard adguard;
    auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(first->Clone(name.c_str())));
    out->SetDirectory(nullptr);
    out->ResetBit(kCanDelete);
    ensure_sumw2(*out);
    out->Reset("ICES");

    for (auto *p : parts)
    {
        if (p)
        {
            out->Add(p);
        }
    }

    if (cfg.fold_overflow)
    {
        const auto [first_bin, last_bin] = interior_bins_for_overflow(*out, cfg.edge_bins);
        // Fold into first/last *interior* bin so edge padding bins can remain empty.
        fold_overflow_into(*out, first_bin, last_bin);
    }
    return out;
}

std::vector<double> AdaptiveBinningService::edges_min_stat(const TH1D &fine,
                                                           const MinStatConfig &cfg) const
{
    if (!cfg.enabled)
    {
        return {};
    }

    if (!(cfg.min_sumw > 0.0) && !(cfg.max_rel_err > 0.0))
    {
        auto edges = uniform_edges_from(fine);
        log_adaptive_bin_sizes(fine.GetName(), edges);
        return edges;
    }

    const TH1D *hptr = &fine;
    std::unique_ptr<TH1D> tmp;
    if (cfg.fold_overflow)
    {
        const std::string tname = std::string(fine.GetName()) + "_minstat_tmp";
        AddDirectoryGuard adguard;
        tmp.reset(static_cast<TH1D *>(fine.Clone(tname.c_str())));
        tmp->SetDirectory(nullptr);
        tmp->ResetBit(kCanDelete);
        ensure_sumw2(*tmp);
        const auto [first_bin, last_bin] = interior_bins_for_overflow(*tmp, cfg.edge_bins);
        fold_overflow_into(*tmp, first_bin, last_bin);
        hptr = tmp.get();
    }
    const TH1D &h = *hptr;

    const int nb = h.GetNbinsX();
    const auto *ax = h.GetXaxis();

    std::vector<double> edges;
    edges.reserve(static_cast<std::size_t>(nb) + 1);

    edges.push_back(ax->GetBinLowEdge(1));

    struct Stat
    {
        double sumw = 0.0;
        double sumw2 = 0.0;
    };
    std::vector<Stat> stats;
    stats.reserve(static_cast<std::size_t>(nb));

    const int edge_bins = std::max(0, cfg.edge_bins);
    const int left_keep = std::min(edge_bins, nb);
    const int right_keep = std::min(edge_bins, nb - left_keep);

    // Keep N fixed-width edge bins on each side (constant small size).
    // These are defined in terms of the *fine* binning.
    for (int i = 1; i <= left_keep; ++i)
    {
        const double up = ax->GetBinUpEdge(i);
        if (up > edges.back() + kEdgeEps)
        {
            edges.push_back(up);
        }
    }

    const int first_bin = left_keep + 1;
    const int last_bin = nb - right_keep;

    double acc_w = 0.0;
    double acc_w2 = 0.0;

    for (int i = first_bin; i <= last_bin; ++i)
    {
        const double w = h.GetBinContent(i);
        acc_w += w;
        // Use explicit sumw2 so this always reflects event weights (and their normalisation),
        // rather than relying on whatever GetBinError currently represents.
        acc_w2 += bin_sumw2(h, i);

        if (pass_bin(acc_w, acc_w2, cfg))
        {
            const double up = ax->GetBinUpEdge(i);
            if (up > edges.back() + kEdgeEps)
            {
                edges.push_back(up);
                stats.push_back(Stat{acc_w, acc_w2});
            }
            acc_w = 0.0;
            acc_w2 = 0.0;
        }
    }

    if (first_bin <= last_bin)
    {
        const double interior_xmax = ax->GetBinUpEdge(last_bin);
        if (edges.back() < interior_xmax - kEdgeEps)
        {
            edges.push_back(interior_xmax);
            stats.push_back(Stat{acc_w, acc_w2});
        }
    }

    // Ensure the last *interior* bin passes; do not merge into right edge padding bins.
    while (stats.size() >= 2)
    {
        const auto &last = stats.back();
        if (pass_bin(last.sumw, last.sumw2, cfg))
        {
            break;
        }

        stats[stats.size() - 2].sumw += last.sumw;
        stats[stats.size() - 2].sumw2 += last.sumw2;

        // Remove the boundary between the last two interior bins:
        // edges currently ends with the interior xmax, so the boundary to remove is edges.end()-2.
        if (edges.size() >= static_cast<std::size_t>(left_keep) + 2u)
        {
            edges.erase(edges.end() - 2);
        }
        stats.pop_back();
    }

    // Append fixed-width right edge bins.
    for (int i = last_bin + 1; i <= nb; ++i)
    {
        const double up = ax->GetBinUpEdge(i);
        if (up > edges.back() + kEdgeEps)
        {
            edges.push_back(up);
        }
    }

    edges.erase(std::unique(edges.begin(), edges.end(),
                            [](double a, double b) { return std::abs(a - b) <= kEdgeEps; }),
                edges.end());

    if (edges.size() < 2)
    {
        edges.clear();
        edges.push_back(ax->GetXmin());
        edges.push_back(ax->GetXmax());
    }

    log_adaptive_bin_sizes(fine.GetName(), edges);
    return edges;
}

std::unique_ptr<TH1D> AdaptiveBinningService::rebin_to_edges(const TH1D &h,
                                                             const std::vector<double> &edges,
                                                             std::string_view new_name,
                                                             const MinStatConfig &cfg) const
{
    const std::string out_name(new_name);

    if (edges.size() < 2)
    {
        auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(h.Clone(out_name.c_str())));
        out->SetDirectory(nullptr);
        ensure_sumw2(*out);
        return out;
    }

    const std::string tmp_name = out_name + "_tmp";
    AddDirectoryGuard adguard;
    auto tmp = std::unique_ptr<TH1D>(static_cast<TH1D *>(h.Clone(tmp_name.c_str())));
    tmp->SetDirectory(nullptr);
    tmp->ResetBit(kCanDelete);
    ensure_sumw2(*tmp);

    if (cfg.fold_overflow)
    {
        const auto [first_bin, last_bin] = interior_bins_for_overflow(*tmp, cfg.edge_bins);
        fold_overflow_into(*tmp, first_bin, last_bin);
    }

    const int nnew = static_cast<int>(edges.size()) - 1;
    TH1 *reb = tmp->Rebin(nnew, out_name.c_str(), edges.data());
    if (!reb)
    {
        // ROOT couldn't rebin (invalid edges, etc.). Return the tmp clone (renamed) so callers
        // still get a valid histogram instead of crashing.
        tmp->SetName(out_name.c_str());
        tmp->SetDirectory(nullptr);
        tmp->ResetBit(kCanDelete);
        ensure_sumw2(*tmp);
        return tmp;
    }

    std::unique_ptr<TH1D> out;
    if (reb == tmp.get())
    {
        // ROOT rebinned in place. Transfer ownership out of tmp to avoid double-free.
        out.reset(tmp.release());
    }
    else
    {
        out.reset(static_cast<TH1D *>(reb));
        // tmp is just a helper clone; let it go now.
        tmp.reset();
    }

    out->SetDirectory(nullptr);
    out->ResetBit(kCanDelete);
    ensure_sumw2(*out);
    return out;
}

std::unique_ptr<TH1D> AdaptiveBinningService::sum_hists(const std::vector<const TH1D *> &parts,
                                                        std::string_view new_name,
                                                        bool fold_overflow) const
{
    MinStatConfig cfg;
    cfg.fold_overflow = fold_overflow;
    cfg.edge_bins = 0;
    return sum_hists(parts, new_name, cfg);
}

std::unique_ptr<TH1D> AdaptiveBinningService::rebin_to_edges(const TH1D &h,
                                                             const std::vector<double> &edges,
                                                             std::string_view new_name,
                                                             bool fold_overflow) const
{
    MinStatConfig cfg;
    cfg.fold_overflow = fold_overflow;
    cfg.edge_bins = 0;
    return rebin_to_edges(h, edges, new_name, cfg);
}

} // namespace nu
