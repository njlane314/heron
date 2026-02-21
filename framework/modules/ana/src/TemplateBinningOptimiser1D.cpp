/* -- C++ -- */
/**
 *  @file  framework/ana/src/TemplateBinningOptimiser1D.cpp
 *
 *  @brief Implementation for template-based 1D binning optimisation.
 */

#include "TemplateBinningOptimizer1D.hh"

#include <TAxis.h>
#include <TDecompSVD.h>
#include <TH1.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace
{

std::vector<double> extract_edges_x(const TH1 &h)
{
    const int n = h.GetNbinsX();
    std::vector<double> edges(n + 1);
    for (int i = 1; i <= n; ++i)
        edges[i - 1] = h.GetXaxis()->GetBinLowEdge(i);
    edges[n] = h.GetXaxis()->GetBinUpEdge(n);
    return edges;
}

bool same_binning_x(const TH1 &a, const TH1 &b, const double eps = 1e-12)
{
    if (a.GetNbinsX() != b.GetNbinsX())
        return false;

    const int n = a.GetNbinsX();
    for (int i = 1; i <= n; ++i)
    {
        const double a_low = a.GetXaxis()->GetBinLowEdge(i);
        const double b_low = b.GetXaxis()->GetBinLowEdge(i);
        if (std::abs(a_low - b_low) > eps)
            return false;
    }

    const double a_high = a.GetXaxis()->GetBinUpEdge(n);
    const double b_high = b.GetXaxis()->GetBinUpEdge(n);
    return (std::abs(a_high - b_high) <= eps);
}

void add_sym_in_place(TMatrixDSym &a, const TMatrixDSym &b, const double scale)
{
    const int n = a.GetNrows();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
            a(i, j) += scale * b(i, j);
    }
}

bool invert_symmetric_svd(const TMatrixDSym &a, TMatrixDSym &a_inverse)
{
    const int n = a.GetNrows();
    TMatrixD dense(n, n);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            dense(i, j) = a(i, j);
    }

    TDecompSVD svd(dense);
    Bool_t ok = kFALSE;
    TMatrixD inverse = svd.Invert(ok);
    if (!ok)
        return false;

    a_inverse.ResizeTo(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            const double sym = 0.5 * (inverse(i, j) + inverse(j, i));
            a_inverse(i, j) = sym;
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
    int n_channel = 0;
    int n_parameter = 0;

    std::vector<double> edges;
    std::vector<std::string> parameter_names;
    std::vector<double> prior_sigmas;
    int poi_index = 0;

    std::vector<double> mu;
    std::vector<double> var;
    std::vector<double> dmu;

    double get_mu(const int c, const int i) const { return mu[c * n_fine + i]; }
    double get_var(const int c, const int i) const { return var[c * n_fine + i]; }
    double get_dmu(const int c, const int p, const int i) const { return dmu[(c * n_parameter + p) * n_fine + i]; }
};

struct BinState
{
    int low = 0;
    int high = 0;

    std::vector<double> mu;
    std::vector<double> var;
    std::vector<double> dmu;

    TMatrixDSym fisher;
};

TMatrixDSym fisher_from_bin_sums(const BinState &bin,
                                 const int n_channel,
                                 const int n_parameter,
                                 const double mu_floor)
{
    TMatrixDSym fisher(n_parameter);

    for (int c = 0; c < n_channel; ++c)
    {
        const double mu = bin.mu[c];
        const double denom = std::max(mu, mu_floor);
        if (denom <= 0.0)
            continue;

        for (int a = 0; a < n_parameter; ++a)
        {
            const double da = bin.dmu[c * n_parameter + a];
            for (int b = 0; b <= a; ++b)
            {
                const double db = bin.dmu[c * n_parameter + b];
                fisher(a, b) += (da * db) / denom;
            }
        }
    }

    return fisher;
}

double sigma_poi_from_fisher(const TMatrixDSym &fisher_without_priors,
                             const std::vector<double> &prior_sigmas,
                             const int poi_index,
                             const bool profile_nuisances)
{
    const int n_parameter = fisher_without_priors.GetNrows();
    if (n_parameter <= 0)
        return std::numeric_limits<double>::infinity();

    if (!profile_nuisances)
    {
        double i_pp = fisher_without_priors(poi_index, poi_index);
        if (prior_sigmas[poi_index] > 0.0)
            i_pp += 1.0 / (prior_sigmas[poi_index] * prior_sigmas[poi_index]);

        if (!(i_pp > 0.0) || !std::isfinite(i_pp))
            return std::numeric_limits<double>::infinity();
        return 1.0 / std::sqrt(i_pp);
    }

    TMatrixDSym total_fisher(fisher_without_priors);
    for (int a = 0; a < n_parameter; ++a)
    {
        if (prior_sigmas[a] > 0.0)
            total_fisher(a, a) += 1.0 / (prior_sigmas[a] * prior_sigmas[a]);
    }

    TMatrixDSym covariance(n_parameter);
    if (!invert_symmetric_svd(total_fisher, covariance))
        return std::numeric_limits<double>::infinity();

    const double variance = covariance(poi_index, poi_index);
    if (!(variance > 0.0) || !std::isfinite(variance))
        return std::numeric_limits<double>::infinity();

    return std::sqrt(variance);
}

ConstraintEval evaluate_constraints(const BinState &bin,
                                    const TemplateBinningOptimizer1D::Config &cfg,
                                    const int n_channel,
                                    const double mu_floor_for_rel)
{
    ConstraintEval out;

    if (cfg.require_per_channel_constraints)
    {
        double worst_penalty = 0.0;
        double mu_sum = 0.0;
        double rel_worst = 0.0;

        for (int c = 0; c < n_channel; ++c)
        {
            const double mu = bin.mu[c];
            const double rel = (mu > 0.0)
                                   ? std::sqrt(std::max(bin.var[c], 0.0)) / std::max(mu, mu_floor_for_rel)
                                   : std::numeric_limits<double>::infinity();

            mu_sum += mu;
            rel_worst = std::max(rel_worst, rel);

            double penalty = 0.0;
            if (mu < cfg.mu_min)
            {
                out.passes = false;
                if (cfg.mu_min > 0.0)
                    penalty += (cfg.mu_min - mu) / cfg.mu_min;
                else
                    penalty += 1.0;
            }
            if (cfg.rel_mc_max > 0.0 && rel > cfg.rel_mc_max)
            {
                out.passes = false;
                penalty += (rel - cfg.rel_mc_max) / cfg.rel_mc_max;
            }

            worst_penalty = std::max(worst_penalty, penalty);
        }

        out.mu_sum = mu_sum;
        out.rel_mc_worst = rel_worst;
        out.penalty = 1000.0 * worst_penalty;
        return out;
    }

    double mu_sum = 0.0;
    double var_sum = 0.0;
    for (int c = 0; c < n_channel; ++c)
    {
        mu_sum += bin.mu[c];
        var_sum += bin.var[c];
    }

    const double rel = (mu_sum > 0.0)
                           ? std::sqrt(std::max(var_sum, 0.0)) / std::max(mu_sum, mu_floor_for_rel)
                           : std::numeric_limits<double>::infinity();

    out.mu_sum = mu_sum;
    out.rel_mc_worst = rel;

    double penalty = 0.0;
    if (mu_sum < cfg.mu_min)
    {
        out.passes = false;
        if (cfg.mu_min > 0.0)
            penalty += (cfg.mu_min - mu_sum) / cfg.mu_min;
        else
            penalty += 1.0;
    }

    if (cfg.rel_mc_max > 0.0 && rel > cfg.rel_mc_max)
    {
        out.passes = false;
        penalty += (rel - cfg.rel_mc_max) / cfg.rel_mc_max;
    }

    out.penalty = 1000.0 * penalty;
    return out;
}

FineCache build_fine_cache(const std::vector<TemplateBinningOptimizer1D::Channel> &channels,
                           const TemplateBinningOptimizer1D::Config &cfg)
{
    if (channels.empty())
        throw std::runtime_error("optimise() called with zero channels");

    const TH1 *p_histogram_0 = channels.front().p_nominal;
    if (!p_histogram_0)
        throw std::runtime_error("channel[0] has null nominal histogram");

    for (size_t c = 0; c < channels.size(); ++c)
    {
        if (!channels[c].p_nominal)
            throw std::runtime_error("a channel has null nominal histogram");

        if (!same_binning_x(*p_histogram_0, *channels[c].p_nominal))
            throw std::runtime_error("channels have different nominal binnings; shared-edge optimisation requires identical fine binning");
    }

    const auto &parameters_0 = channels[0].parameters;
    if (parameters_0.empty())
        throw std::runtime_error("channel[0] has no parameters; need at least a POI for objective");

    for (size_t c = 1; c < channels.size(); ++c)
    {
        const auto &parameters = channels[c].parameters;
        if (parameters.size() != parameters_0.size())
            throw std::runtime_error("channels have different parameter counts");

        for (size_t p = 0; p < parameters_0.size(); ++p)
        {
            if (parameters[p].name != parameters_0[p].name)
                throw std::runtime_error("channels have different parameter name/order");
            if (parameters[p].is_poi != parameters_0[p].is_poi)
                throw std::runtime_error("channels have inconsistent is_poi flags");
            if (std::abs(parameters[p].prior_sigma - parameters_0[p].prior_sigma) > 0.0)
                throw std::runtime_error("channels have inconsistent prior_sigma for a parameter");
        }
    }

    FineCache cache;
    cache.n_fine = p_histogram_0->GetNbinsX();
    cache.n_channel = static_cast<int>(channels.size());
    cache.n_parameter = static_cast<int>(parameters_0.size());
    cache.edges = extract_edges_x(*p_histogram_0);

    cache.parameter_names.resize(cache.n_parameter);
    cache.prior_sigmas.resize(cache.n_parameter, 0.0);

    int poi_index = -1;
    for (int p = 0; p < cache.n_parameter; ++p)
    {
        cache.parameter_names[p] = parameters_0[p].name;
        cache.prior_sigmas[p] = parameters_0[p].prior_sigma;
        if (parameters_0[p].is_poi)
        {
            if (poi_index >= 0)
                throw std::runtime_error("multiple parameters marked is_poi=true; require exactly one");
            poi_index = p;
        }
    }
    if (poi_index < 0)
        poi_index = 0;
    cache.poi_index = poi_index;

    for (int c = 0; c < cache.n_channel; ++c)
    {
        const TH1 &nominal = *channels[c].p_nominal;

        for (int p = 0; p < cache.n_parameter; ++p)
        {
            const auto &parameter = channels[c].parameters[p];
            if (parameter.p_derivative)
            {
                if (!same_binning_x(nominal, *parameter.p_derivative))
                    throw std::runtime_error("derivative histogram has different binning than nominal");
            }
            else
            {
                if (!parameter.p_up || !parameter.p_down)
                    throw std::runtime_error("parameter must provide either derivative or (up and down)");
                if (!same_binning_x(nominal, *parameter.p_up) || !same_binning_x(nominal, *parameter.p_down))
                    throw std::runtime_error("up/down histogram has different binning than nominal");
                if (!(parameter.step > 0.0))
                    throw std::runtime_error("parameter step must be > 0 for up/down finite difference");
            }
        }
    }

    cache.mu.assign(cache.n_channel * cache.n_fine, 0.0);
    cache.var.assign(cache.n_channel * cache.n_fine, 0.0);
    cache.dmu.assign(cache.n_channel * cache.n_parameter * cache.n_fine, 0.0);

    for (int c = 0; c < cache.n_channel; ++c)
    {
        const TH1 &nominal = *channels[c].p_nominal;

        for (int i = 1; i <= cache.n_fine; ++i)
        {
            const double mu = nominal.GetBinContent(i);
            const double err = nominal.GetBinError(i);
            cache.mu[c * cache.n_fine + (i - 1)] = mu;
            cache.var[c * cache.n_fine + (i - 1)] = err * err;

            for (int p = 0; p < cache.n_parameter; ++p)
            {
                const auto &parameter = channels[c].parameters[p];
                double derivative = 0.0;

                if (parameter.p_derivative)
                    derivative = parameter.p_derivative->GetBinContent(i);
                else
                {
                    const double up = parameter.p_up->GetBinContent(i);
                    const double down = parameter.p_down->GetBinContent(i);
                    derivative = (up - down) / (2.0 * parameter.step);
                }

                cache.dmu[(c * cache.n_parameter + p) * cache.n_fine + (i - 1)] = derivative;
            }
        }
    }

    if (cfg.verbose && cfg.p_log)
    {
        (*cfg.p_log) << "[TemplateBinningOptimizer1D] fine bins=" << cache.n_fine
                     << " channels=" << cache.n_channel
                     << " parameters=" << cache.n_parameter
                     << " POI=" << cache.parameter_names[cache.poi_index] << "\n";
    }

    return cache;
}

BinState make_fine_bin_state(const FineCache &cache,
                             const int i_fine,
                             const double mu_floor_for_objective)
{
    BinState state;
    state.low = i_fine;
    state.high = i_fine;

    state.mu.assign(cache.n_channel, 0.0);
    state.var.assign(cache.n_channel, 0.0);
    state.dmu.assign(cache.n_channel * cache.n_parameter, 0.0);

    for (int c = 0; c < cache.n_channel; ++c)
    {
        state.mu[c] = cache.get_mu(c, i_fine);
        state.var[c] = cache.get_var(c, i_fine);
        for (int p = 0; p < cache.n_parameter; ++p)
            state.dmu[c * cache.n_parameter + p] = cache.get_dmu(c, p, i_fine);
    }

    state.fisher = fisher_from_bin_sums(state, cache.n_channel, cache.n_parameter, mu_floor_for_objective);
    return state;
}

BinState merge_bins(const BinState &left,
                    const BinState &right,
                    const int n_channel,
                    const int n_parameter,
                    const double mu_floor_for_objective)
{
    BinState merged;
    merged.low = std::min(left.low, right.low);
    merged.high = std::max(left.high, right.high);

    merged.mu.resize(n_channel);
    merged.var.resize(n_channel);
    merged.dmu.resize(n_channel * n_parameter);

    for (int c = 0; c < n_channel; ++c)
    {
        merged.mu[c] = left.mu[c] + right.mu[c];
        merged.var[c] = left.var[c] + right.var[c];
        for (int p = 0; p < n_parameter; ++p)
            merged.dmu[c * n_parameter + p] = left.dmu[c * n_parameter + p] + right.dmu[c * n_parameter + p];
    }

    merged.fisher = fisher_from_bin_sums(merged, n_channel, n_parameter, mu_floor_for_objective);
    return merged;
}

} // namespace

TemplateBinningOptimizer1D::TemplateBinningOptimizer1D(Config cfg) : m_cfg(std::move(cfg)) {}

TemplateBinningOptimizer1D::Result TemplateBinningOptimizer1D::optimise(const Channel &channel) const
{
    return optimise(std::vector<Channel>{channel});
}

TemplateBinningOptimizer1D::Result TemplateBinningOptimizer1D::optimise(const std::vector<Channel> &channels) const
{
    FineCache cache = build_fine_cache(channels, m_cfg);

    std::vector<BinState> bins;
    bins.reserve(cache.n_fine);
    for (int i = 0; i < cache.n_fine; ++i)
        bins.push_back(make_fine_bin_state(cache, i, m_cfg.mu_floor_for_objective));

    TMatrixDSym total_fisher(cache.n_parameter);
    for (const auto &bin : bins)
        add_sym_in_place(total_fisher, bin.fisher, +1.0);

    double sigma_current = sigma_poi_from_fisher(total_fisher, cache.prior_sigmas, cache.poi_index, m_cfg.profile_nuisances);

    auto need_more_merging = [&](const std::vector<BinState> &states) {
        bool any_failing = false;
        for (const auto &state : states)
        {
            const auto eval = evaluate_constraints(state, m_cfg, cache.n_channel, m_cfg.mu_floor_for_objective);
            if (!eval.passes)
            {
                any_failing = true;
                break;
            }
        }

        const bool too_many_bins = (m_cfg.max_bins > 0) && (static_cast<int>(states.size()) > m_cfg.max_bins);
        return any_failing || too_many_bins;
    };

    int iteration = 0;
    while (bins.size() > 1 && need_more_merging(bins))
    {
        ++iteration;

        std::vector<bool> fails(bins.size(), false);
        bool any_failing = false;
        for (size_t i = 0; i < bins.size(); ++i)
        {
            const auto eval = evaluate_constraints(bins[i], m_cfg, cache.n_channel, m_cfg.mu_floor_for_objective);
            fails[i] = !eval.passes;
            any_failing = any_failing || fails[i];
        }

        int best_index = -1;
        double best_cost = std::numeric_limits<double>::infinity();
        double best_sigma = sigma_current;
        BinState best_merged;

        for (int k = 0; k < static_cast<int>(bins.size()) - 1; ++k)
        {
            const bool touches_failing_bin = fails[k] || fails[k + 1];
            if (any_failing && !touches_failing_bin)
                continue;

            BinState merged = merge_bins(bins[k], bins[k + 1], cache.n_channel, cache.n_parameter, m_cfg.mu_floor_for_objective);

            TMatrixDSym candidate_fisher(total_fisher);
            add_sym_in_place(candidate_fisher, bins[k].fisher, -1.0);
            add_sym_in_place(candidate_fisher, bins[k + 1].fisher, -1.0);
            add_sym_in_place(candidate_fisher, merged.fisher, +1.0);

            const double sigma_candidate =
                sigma_poi_from_fisher(candidate_fisher, cache.prior_sigmas, cache.poi_index, m_cfg.profile_nuisances);

            const auto merged_eval = evaluate_constraints(merged, m_cfg, cache.n_channel, m_cfg.mu_floor_for_objective);
            const double width = cache.edges[merged.high + 1] - cache.edges[merged.low];

            const double cost =
                (sigma_candidate - sigma_current) +
                merged_eval.penalty +
                m_cfg.width_penalty * width;

            if (cost < best_cost)
            {
                best_cost = cost;
                best_index = k;
                best_sigma = sigma_candidate;
                best_merged = std::move(merged);
            }
        }

        if (best_index < 0)
        {
            for (int k = 0; k < static_cast<int>(bins.size()) - 1; ++k)
            {
                BinState merged = merge_bins(bins[k], bins[k + 1], cache.n_channel, cache.n_parameter, m_cfg.mu_floor_for_objective);

                TMatrixDSym candidate_fisher(total_fisher);
                add_sym_in_place(candidate_fisher, bins[k].fisher, -1.0);
                add_sym_in_place(candidate_fisher, bins[k + 1].fisher, -1.0);
                add_sym_in_place(candidate_fisher, merged.fisher, +1.0);

                const double sigma_candidate =
                    sigma_poi_from_fisher(candidate_fisher, cache.prior_sigmas, cache.poi_index, m_cfg.profile_nuisances);
                const auto merged_eval = evaluate_constraints(merged, m_cfg, cache.n_channel, m_cfg.mu_floor_for_objective);
                const double width = cache.edges[merged.high + 1] - cache.edges[merged.low];

                const double cost =
                    (sigma_candidate - sigma_current) +
                    merged_eval.penalty +
                    m_cfg.width_penalty * width;

                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_index = k;
                    best_sigma = sigma_candidate;
                    best_merged = std::move(merged);
                }
            }
        }

        if (best_index < 0)
            break;

        add_sym_in_place(total_fisher, bins[best_index].fisher, -1.0);
        add_sym_in_place(total_fisher, bins[best_index + 1].fisher, -1.0);
        add_sym_in_place(total_fisher, best_merged.fisher, +1.0);

        bins[best_index] = std::move(best_merged);
        bins.erase(bins.begin() + best_index + 1);

        sigma_current = best_sigma;

        if (m_cfg.verbose && m_cfg.p_log)
        {
            (*m_cfg.p_log) << "[TemplateBinningOptimizer1D] iter=" << iteration
                           << " bins=" << bins.size()
                           << " expected_sigma_poi=" << sigma_current
                           << "\n";
        }
    }

    Result out;
    out.edges.reserve(bins.size() + 1);
    out.edges.push_back(cache.edges[bins.front().low]);
    for (const auto &bin : bins)
        out.edges.push_back(cache.edges[bin.high + 1]);

    out.expected_sigma_poi = sigma_current;

    out.bins.reserve(bins.size());
    for (const auto &bin : bins)
    {
        BinReport report;
        report.low = cache.edges[bin.low];
        report.high = cache.edges[bin.high + 1];

        const auto eval = evaluate_constraints(bin, m_cfg, cache.n_channel, m_cfg.mu_floor_for_objective);
        report.mu_sum = eval.mu_sum;
        report.rel_mc_worst = eval.rel_mc_worst;
        report.passes_constraints = eval.passes;

        out.bins.push_back(report);
    }

    return out;
}
