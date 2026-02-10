/* -- C++ -- */
/**
 *  @file  plot/src/AdaptiveBinningService.cc
 *
 *  @brief Adaptive binning service implementation.
 */

#include "AdaptiveBinningService.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

#include <TAxis.h>
#include <TH1D.h>

#include "PlotDescriptors.hh"

namespace nu
{
namespace
{
constexpr double kEdgeEps = 1e-12;

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
    {
        return;
    }

    std::ostringstream msg;
    msg << "[AdaptiveBinningService] Adaptive bins settled for '" << hist_name << "': "
        << (edges.size() - 1) << " bins; widths [";

    for (std::size_t i = 1; i < edges.size(); ++i)
    {
        if (i > 1)
        {
            msg << ", ";
        }
        msg << (edges[i] - edges[i - 1]);
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

    {
        const double c0 = h.GetBinContent(0);
        const double e0 = h.GetBinError(0);
        const double c1 = h.GetBinContent(1);
        const double e1 = h.GetBinError(1);
        h.SetBinContent(1, c1 + c0);
        h.SetBinError(1, std::hypot(e1, e0));
        h.SetBinContent(0, 0.0);
        h.SetBinError(0, 0.0);
    }

    {
        const double co = h.GetBinContent(nb + 1);
        const double eo = h.GetBinError(nb + 1);
        const double cn = h.GetBinContent(nb);
        const double en = h.GetBinError(nb);
        h.SetBinContent(nb, cn + co);
        h.SetBinError(nb, std::hypot(en, eo));
        h.SetBinContent(nb + 1, 0.0);
        h.SetBinError(nb + 1, 0.0);
    }
}

std::unique_ptr<TH1D> AdaptiveBinningService::sum_hists(const std::vector<const TH1D *> &parts,
                                                        std::string_view new_name,
                                                        bool do_fold_overflow) const
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
    auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(first->Clone(name.c_str())));
    out->SetDirectory(nullptr);
    out->Sumw2(kTRUE);
    out->Reset("ICES");

    for (auto *p : parts)
    {
        if (p)
        {
            out->Add(p);
        }
    }

    if (do_fold_overflow)
    {
        fold_overflow(*out);
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
        tmp.reset(static_cast<TH1D *>(fine.Clone(tname.c_str())));
        tmp->SetDirectory(nullptr);
        tmp->Sumw2(kTRUE);
        fold_overflow(*tmp);
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

    double acc_w = 0.0;
    double acc_w2 = 0.0;

    for (int i = 1; i <= nb; ++i)
    {
        const double w = h.GetBinContent(i);
        const double e = h.GetBinError(i);
        acc_w += w;
        acc_w2 += e * e;

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

    const double xmax = ax->GetBinUpEdge(nb);
    if (edges.back() < xmax - kEdgeEps)
    {
        edges.push_back(xmax);
        stats.push_back(Stat{acc_w, acc_w2});
    }

    while (stats.size() >= 2)
    {
        const auto &last = stats.back();
        if (pass_bin(last.sumw, last.sumw2, cfg))
        {
            break;
        }

        stats[stats.size() - 2].sumw += last.sumw;
        stats[stats.size() - 2].sumw2 += last.sumw2;
        stats.pop_back();

        if (edges.size() >= 3)
        {
            edges.erase(edges.end() - 2);
        }
        else
        {
            break;
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
                                                             bool do_fold_overflow) const
{
    const std::string out_name(new_name);

    if (edges.size() < 2)
    {
        auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(h.Clone(out_name.c_str())));
        out->SetDirectory(nullptr);
        out->Sumw2(kTRUE);
        return out;
    }

    const std::string tmp_name = out_name + "_tmp";
    auto tmp = std::unique_ptr<TH1D>(static_cast<TH1D *>(h.Clone(tmp_name.c_str())));
    tmp->SetDirectory(nullptr);
    tmp->Sumw2(kTRUE);

    if (do_fold_overflow)
    {
        fold_overflow(*tmp);
    }

    const int nnew = static_cast<int>(edges.size()) - 1;
    auto *reb = tmp->Rebin(nnew, out_name.c_str(), edges.data());

    auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(reb));
    out->SetDirectory(nullptr);
    out->Sumw2(kTRUE);
    return out;
}

} // namespace nu
