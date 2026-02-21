#include "TemplateBinningOptimiser1D.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <TAxis.h>
#include <TDecompSVD.h>
#include <TMatrixD.h>

namespace nu
{
namespace binning
{
namespace
{

std::vector<double> extract_edges_x(const TH1 &h)
{
    const int n = h.GetNbinsX();
    std::vector<double> edges(n + 1);
    for (int i = 1; i <= n; ++i)
    {
        edges[i - 1] = h.GetXaxis()->GetBinLowEdge(i);
    }
    edges[n] = h.GetXaxis()->GetBinUpEdge(n);
    return edges;
}

bool same_binning_x(const TH1 &a, const TH1 &b, double eps = 1e-12)
{
    if (a.GetNbinsX() != b.GetNbinsX()) return false;

    const int n = a.GetNbinsX();
    for (int i = 1; i <= n; ++i)
    {
        const double al = a.GetXaxis()->GetBinLowEdge(i);
        const double bl = b.GetXaxis()->GetBinLowEdge(i);
        if (std::abs(al - bl) > eps) return false;
    }

    const double ah = a.GetXaxis()->GetBinUpEdge(n);
    const double bh = b.GetXaxis()->GetBinUpEdge(n);
    return std::abs(ah - bh) <= eps;
}

void add_sym_in_place(TMatrixDSym &a, const TMatrixDSym &b, double scale)
{
    const int n = a.GetNrows();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            a(i, j) += scale * b(i, j);
        }
    }
}

bool invert_symmetric_svd(const TMatrixDSym &a, TMatrixDSym &inv)
{
    const int n = a.GetNrows();
    TMatrixD ad(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            ad(i, j) = a(i, j);
        }
    }

    TDecompSVD svd(ad);
    Bool_t ok = kFALSE;
    TMatrixD ad_inv = svd.Invert(ok);
    if (!ok) return false;

    inv.ResizeTo(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            inv(i, j) = 0.5 * (ad_inv(i, j) + ad_inv(j, i));
        }
    }
    return true;
}

struct ConstraintEval
{
    bool passes = true;
    double penalty = 0.0;
    double mu_sum = 0.0;
    double rel_mc_worst = 0.0;
};

struct FineCache
{
    int n_fine = 0;
    int n_chan = 0;
    int n_par = 0;

    std::vector<double> edges;
    std::vector<std::string> par_names;
    std::vector<double> prior_sigma;
    int poi_index = 0;

    std::vector<double> mu;
    std::vector<double> var;
    std::vector<double> dmu;

    double Mu(int c, int i) const { return mu[c * n_fine + i]; }
    double Var(int c, int i) const { return var[c * n_fine + i]; }
    double Dmu(int c, int a, int i) const { return dmu[(c * n_par + a) * n_fine + i]; }
};

struct BinState
{
    int lo = 0;
    int hi = 0;

    std::vector<double> mu;
    std::vector<double> var;
    std::vector<double> dmu;

    TMatrixDSym fisher;
};

TMatrixDSym fisher_from_bin_sums(const BinState &b, int n_chan, int n_par, double mu_floor_for_objective)
{
    TMatrixDSym info(n_par);
    for (int c = 0; c < n_chan; ++c)
    {
        const double mu = b.mu[c];
        const double denom = std::max(mu, mu_floor_for_objective);
        if (denom <= 0.0) continue;

        for (int a = 0; a < n_par; ++a)
        {
            const double da = b.dmu[c * n_par + a];
            for (int bb = 0; bb <= a; ++bb)
            {
                const double db = b.dmu[c * n_par + bb];
                info(a, bb) += (da * db) / denom;
            }
        }
    }
    return info;
}

double sigma_poi_from_fisher(const TMatrixDSym &fisher_no_priors,
                             const std::vector<double> &prior_sigma,
                             int poi_index,
                             bool profile_nuisances)
{
    const int p = fisher_no_priors.GetNrows();
    if (p <= 0) return std::numeric_limits<double>::infinity();

    if (!profile_nuisances)
    {
        double poi_info = fisher_no_priors(poi_index, poi_index);
        if (prior_sigma[poi_index] > 0.0)
        {
            poi_info += 1.0 / (prior_sigma[poi_index] * prior_sigma[poi_index]);
        }
        if (!(poi_info > 0.0) || !std::isfinite(poi_info)) return std::numeric_limits<double>::infinity();
        return 1.0 / std::sqrt(poi_info);
    }

    TMatrixDSym total(fisher_no_priors);
    for (int a = 0; a < p; ++a)
    {
        if (prior_sigma[a] > 0.0)
        {
            total(a, a) += 1.0 / (prior_sigma[a] * prior_sigma[a]);
        }
    }

    TMatrixDSym cov(p);
    if (!invert_symmetric_svd(total, cov)) return std::numeric_limits<double>::infinity();

    const double value = cov(poi_index, poi_index);
    if (!(value > 0.0) || !std::isfinite(value)) return std::numeric_limits<double>::infinity();
    return std::sqrt(value);
}

ConstraintEval eval_constraints(const BinState &b,
                                const TemplateBinningOptimiser1D::Config &cfg,
                                int n_chan,
                                double mu_floor_for_rel)
{
    ConstraintEval out;

    if (cfg.require_per_channel_constraints)
    {
        double worst_penalty = 0.0;
        double mu_sum = 0.0;
        double rel_worst = 0.0;

        for (int c = 0; c < n_chan; ++c)
        {
            const double mu = b.mu[c];
            const double rel = (mu > 0.0)
                                 ? std::sqrt(std::max(b.var[c], 0.0)) / std::max(mu, mu_floor_for_rel)
                                 : std::numeric_limits<double>::infinity();

            mu_sum += mu;
            rel_worst = std::max(rel_worst, rel);

            double pen = 0.0;
            if (mu < cfg.mu_min)
            {
                out.passes = false;
                pen += (cfg.mu_min > 0.0) ? (cfg.mu_min - mu) / cfg.mu_min : 1.0;
            }
            if (cfg.rel_mc_max > 0.0 && rel > cfg.rel_mc_max)
            {
                out.passes = false;
                pen += (rel - cfg.rel_mc_max) / cfg.rel_mc_max;
            }

            worst_penalty = std::max(worst_penalty, pen);
        }

        out.mu_sum = mu_sum;
        out.rel_mc_worst = rel_worst;
        out.penalty = 1000.0 * worst_penalty;
        return out;
    }

    double mu_sum = 0.0;
    double var_sum = 0.0;
    for (int c = 0; c < n_chan; ++c)
    {
        mu_sum += b.mu[c];
        var_sum += b.var[c];
    }

    const double rel = (mu_sum > 0.0)
                           ? std::sqrt(std::max(var_sum, 0.0)) / std::max(mu_sum, mu_floor_for_rel)
                           : std::numeric_limits<double>::infinity();

    out.mu_sum = mu_sum;
    out.rel_mc_worst = rel;

    double pen = 0.0;
    if (mu_sum < cfg.mu_min)
    {
        out.passes = false;
        pen += (cfg.mu_min > 0.0) ? (cfg.mu_min - mu_sum) / cfg.mu_min : 1.0;
    }
    if (cfg.rel_mc_max > 0.0 && rel > cfg.rel_mc_max)
    {
        out.passes = false;
        pen += (rel - cfg.rel_mc_max) / cfg.rel_mc_max;
    }

    out.penalty = 1000.0 * pen;
    return out;
}

FineCache build_fine_cache(const std::vector<TemplateBinningOptimiser1D::Channel> &channels,
                           const TemplateBinningOptimiser1D::Config &cfg)
{
    if (channels.empty()) throw std::runtime_error("optimise() called with zero channels");

    const TH1 *h0 = channels[0].nominal;
    if (!h0) throw std::runtime_error("Channel[0] has null nominal histogram");

    for (size_t c = 0; c < channels.size(); ++c)
    {
        if (!channels[c].nominal) throw std::runtime_error("A channel has null nominal histogram");
        if (!same_binning_x(*h0, *channels[c].nominal))
        {
            throw std::runtime_error("Channels have different nominal binnings.");
        }
    }

    const auto &pars0 = channels[0].parameters;
    if (pars0.empty()) throw std::runtime_error("Channel[0] has no parameters.");

    for (size_t c = 1; c < channels.size(); ++c)
    {
        const auto &pc = channels[c].parameters;
        if (pc.size() != pars0.size()) throw std::runtime_error("Channels have different parameter counts.");

        for (size_t a = 0; a < pars0.size(); ++a)
        {
            if (pc[a].name != pars0[a].name) throw std::runtime_error("Channels have different parameter name/order.");
            if (pc[a].is_poi != pars0[a].is_poi) throw std::runtime_error("Channels have inconsistent is_poi flags.");
            if (std::abs(pc[a].prior_sigma - pars0[a].prior_sigma) > 0.0)
            {
                throw std::runtime_error("Channels have inconsistent prior_sigma for a parameter.");
            }
        }
    }

    FineCache cache;
    cache.n_fine = h0->GetNbinsX();
    cache.n_chan = static_cast<int>(channels.size());
    cache.n_par = static_cast<int>(pars0.size());
    cache.edges = extract_edges_x(*h0);

    cache.par_names.resize(cache.n_par);
    cache.prior_sigma.resize(cache.n_par, 0.0);

    int poi_idx = -1;
    for (int a = 0; a < cache.n_par; ++a)
    {
        cache.par_names[a] = pars0[a].name;
        cache.prior_sigma[a] = pars0[a].prior_sigma;
        if (pars0[a].is_poi)
        {
            if (poi_idx >= 0) throw std::runtime_error("Multiple parameters are marked is_poi=true.");
            poi_idx = a;
        }
    }
    if (poi_idx < 0) poi_idx = 0;
    cache.poi_index = poi_idx;

    for (int c = 0; c < cache.n_chan; ++c)
    {
        const TH1 &h_nom = *channels[c].nominal;
        for (int a = 0; a < cache.n_par; ++a)
        {
            const auto &p = channels[c].parameters[a];
            if (p.derivative)
            {
                if (!same_binning_x(h_nom, *p.derivative))
                {
                    throw std::runtime_error("Derivative histogram has different binning than nominal.");
                }
            }
            else
            {
                if (!p.up || !p.down)
                {
                    throw std::runtime_error("Parameter must provide either derivative or both up/down histograms.");
                }
                if (!same_binning_x(h_nom, *p.up) || !same_binning_x(h_nom, *p.down))
                {
                    throw std::runtime_error("Up/down histogram has different binning than nominal.");
                }
                if (!(p.step > 0.0))
                {
                    throw std::runtime_error("Parameter step must be > 0 for up/down finite difference.");
                }
            }
        }
    }

    cache.mu.assign(cache.n_chan * cache.n_fine, 0.0);
    cache.var.assign(cache.n_chan * cache.n_fine, 0.0);
    cache.dmu.assign(cache.n_chan * cache.n_par * cache.n_fine, 0.0);

    for (int c = 0; c < cache.n_chan; ++c)
    {
        const TH1 &h_nom = *channels[c].nominal;
        for (int i = 1; i <= cache.n_fine; ++i)
        {
            const double mu = h_nom.GetBinContent(i);
            const double err = h_nom.GetBinError(i);
            cache.mu[c * cache.n_fine + (i - 1)] = mu;
            cache.var[c * cache.n_fine + (i - 1)] = err * err;

            for (int a = 0; a < cache.n_par; ++a)
            {
                const auto &p = channels[c].parameters[a];
                double d = 0.0;
                if (p.derivative)
                {
                    d = p.derivative->GetBinContent(i);
                }
                else
                {
                    const double up = p.up->GetBinContent(i);
                    const double down = p.down->GetBinContent(i);
                    d = (up - down) / (2.0 * p.step);
                }
                cache.dmu[(c * cache.n_par + a) * cache.n_fine + (i - 1)] = d;
            }
        }
    }

    if (cfg.verbose && cfg.log)
    {
        (*cfg.log) << "[TemplateBinningOptimiser1D] fine bins=" << cache.n_fine
                   << " channels=" << cache.n_chan
                   << " parameters=" << cache.n_par
                   << " POI=" << cache.par_names[cache.poi_index] << "\n";
    }

    return cache;
}

BinState make_fine_bin_state(const FineCache &cache, int i_fine, double mu_floor_for_objective)
{
    BinState b;
    b.lo = i_fine;
    b.hi = i_fine;

    b.mu.assign(cache.n_chan, 0.0);
    b.var.assign(cache.n_chan, 0.0);
    b.dmu.assign(cache.n_chan * cache.n_par, 0.0);

    for (int c = 0; c < cache.n_chan; ++c)
    {
        b.mu[c] = cache.Mu(c, i_fine);
        b.var[c] = cache.Var(c, i_fine);
        for (int a = 0; a < cache.n_par; ++a)
        {
            b.dmu[c * cache.n_par + a] = cache.Dmu(c, a, i_fine);
        }
    }

    b.fisher = fisher_from_bin_sums(b, cache.n_chan, cache.n_par, mu_floor_for_objective);
    return b;
}

BinState merge_bins(const BinState &left,
                    const BinState &right,
                    int n_chan,
                    int n_par,
                    double mu_floor_for_objective)
{
    BinState merged;
    merged.lo = std::min(left.lo, right.lo);
    merged.hi = std::max(left.hi, right.hi);

    merged.mu.resize(n_chan);
    merged.var.resize(n_chan);
    merged.dmu.resize(n_chan * n_par);

    for (int c = 0; c < n_chan; ++c)
    {
        merged.mu[c] = left.mu[c] + right.mu[c];
        merged.var[c] = left.var[c] + right.var[c];
        for (int a = 0; a < n_par; ++a)
        {
            merged.dmu[c * n_par + a] = left.dmu[c * n_par + a] + right.dmu[c * n_par + a];
        }
    }

    merged.fisher = fisher_from_bin_sums(merged, n_chan, n_par, mu_floor_for_objective);
    return merged;
}

} // namespace

TemplateBinningOptimiser1D::TemplateBinningOptimiser1D(Config cfg) : cfg_(std::move(cfg)) {}

TemplateBinningOptimiser1D::Result TemplateBinningOptimiser1D::optimise(const Channel &channel) const
{
    return optimise(std::vector<Channel>{channel});
}

TemplateBinningOptimiser1D::Result TemplateBinningOptimiser1D::optimise(const std::vector<Channel> &channels) const
{
    FineCache cache = build_fine_cache(channels, cfg_);

    std::vector<BinState> bins;
    bins.reserve(cache.n_fine);
    for (int i = 0; i < cache.n_fine; ++i)
    {
        bins.push_back(make_fine_bin_state(cache, i, cfg_.mu_floor_for_objective));
    }

    TMatrixDSym total_fisher(cache.n_par);
    for (const auto &b : bins)
    {
        add_sym_in_place(total_fisher, b.fisher, +1.0);
    }

    double sigma_current = sigma_poi_from_fisher(total_fisher,
                                                 cache.prior_sigma,
                                                 cache.poi_index,
                                                 cfg_.profile_nuisances);

    auto need_more_merging = [&](const std::vector<BinState> &states)
    {
        bool any_fail = false;
        for (const auto &b : states)
        {
            auto eval = eval_constraints(b, cfg_, cache.n_chan, cfg_.mu_floor_for_objective);
            if (!eval.passes)
            {
                any_fail = true;
                break;
            }
        }

        const bool too_many_bins = (cfg_.max_bins > 0) && (static_cast<int>(states.size()) > cfg_.max_bins);
        return any_fail || too_many_bins;
    };

    int iter = 0;
    while (bins.size() > 1 && need_more_merging(bins))
    {
        ++iter;

        std::vector<bool> fails(bins.size(), false);
        bool any_fail = false;
        for (size_t i = 0; i < bins.size(); ++i)
        {
            auto eval = eval_constraints(bins[i], cfg_, cache.n_chan, cfg_.mu_floor_for_objective);
            fails[i] = !eval.passes;
            any_fail = any_fail || fails[i];
        }

        int best_k = -1;
        double best_cost = std::numeric_limits<double>::infinity();
        double best_sigma = sigma_current;
        BinState best_merged;

        for (int k = 0; k < static_cast<int>(bins.size()) - 1; ++k)
        {
            const bool touches_fail = fails[k] || fails[k + 1];
            if (any_fail && !touches_fail) continue;

            BinState merged = merge_bins(bins[k], bins[k + 1], cache.n_chan, cache.n_par, cfg_.mu_floor_for_objective);

            TMatrixDSym candidate_fisher(total_fisher);
            add_sym_in_place(candidate_fisher, bins[k].fisher, -1.0);
            add_sym_in_place(candidate_fisher, bins[k + 1].fisher, -1.0);
            add_sym_in_place(candidate_fisher, merged.fisher, +1.0);

            const double sigma_candidate = sigma_poi_from_fisher(candidate_fisher,
                                                                 cache.prior_sigma,
                                                                 cache.poi_index,
                                                                 cfg_.profile_nuisances);

            auto merged_eval = eval_constraints(merged, cfg_, cache.n_chan, cfg_.mu_floor_for_objective);
            const double width = cache.edges[merged.hi + 1] - cache.edges[merged.lo];

            const double cost = (sigma_candidate - sigma_current) +
                                merged_eval.penalty +
                                cfg_.width_penalty * width;

            if (cost < best_cost)
            {
                best_cost = cost;
                best_k = k;
                best_sigma = sigma_candidate;
                best_merged = std::move(merged);
            }
        }

        if (best_k < 0)
        {
            for (int k = 0; k < static_cast<int>(bins.size()) - 1; ++k)
            {
                BinState merged = merge_bins(bins[k], bins[k + 1], cache.n_chan, cache.n_par, cfg_.mu_floor_for_objective);

                TMatrixDSym candidate_fisher(total_fisher);
                add_sym_in_place(candidate_fisher, bins[k].fisher, -1.0);
                add_sym_in_place(candidate_fisher, bins[k + 1].fisher, -1.0);
                add_sym_in_place(candidate_fisher, merged.fisher, +1.0);

                const double sigma_candidate = sigma_poi_from_fisher(candidate_fisher,
                                                                     cache.prior_sigma,
                                                                     cache.poi_index,
                                                                     cfg_.profile_nuisances);

                auto merged_eval = eval_constraints(merged, cfg_, cache.n_chan, cfg_.mu_floor_for_objective);
                const double width = cache.edges[merged.hi + 1] - cache.edges[merged.lo];

                const double cost = (sigma_candidate - sigma_current) +
                                    merged_eval.penalty +
                                    cfg_.width_penalty * width;

                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_k = k;
                    best_sigma = sigma_candidate;
                    best_merged = std::move(merged);
                }
            }
        }

        if (best_k < 0) break;

        add_sym_in_place(total_fisher, bins[best_k].fisher, -1.0);
        add_sym_in_place(total_fisher, bins[best_k + 1].fisher, -1.0);
        add_sym_in_place(total_fisher, best_merged.fisher, +1.0);

        bins[best_k] = std::move(best_merged);
        bins.erase(bins.begin() + best_k + 1);

        sigma_current = best_sigma;

        if (cfg_.verbose && cfg_.log)
        {
            (*cfg_.log) << "[TemplateBinningOptimiser1D] iter=" << iter
                        << " bins=" << bins.size()
                        << " expected_sigma_poi=" << sigma_current
                        << "\n";
        }
    }

    Result out;
    out.edges.reserve(bins.size() + 1);
    out.edges.push_back(cache.edges[bins.front().lo]);
    for (const auto &b : bins)
    {
        out.edges.push_back(cache.edges[b.hi + 1]);
    }
    out.expected_sigma_poi = sigma_current;

    out.bins.reserve(bins.size());
    for (const auto &b : bins)
    {
        BinReport report;
        report.low = cache.edges[b.lo];
        report.high = cache.edges[b.hi + 1];

        auto eval = eval_constraints(b, cfg_, cache.n_chan, cfg_.mu_floor_for_objective);
        report.mu_sum = eval.mu_sum;
        report.rel_mc_worst = eval.rel_mc_worst;
        report.passes_constraints = eval.passes;

        out.bins.push_back(report);
    }

    return out;
}

} // namespace binning
} // namespace nu
