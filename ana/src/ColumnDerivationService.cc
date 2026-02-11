/* -- C++ -- */
/**
 *  @file  ana/src/ColumnDerivationService.cc
 *
 *  @brief Variable definitions for analysis RDataFrame processing.
 */

#include "ColumnDerivationService.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <ROOT/RVec.hxx>

#include "SelectionService.hh"

namespace
{

ROOT::RVec<float> derive_track_mipness_median_plane_score(
    const ROOT::RVec<float> &t70_u,
    const ROOT::RVec<float> &q50_u,
    const ROOT::RVec<float> &q90_u,
    const ROOT::RVec<float> &t70_v,
    const ROOT::RVec<float> &q50_v,
    const ROOT::RVec<float> &q90_v,
    const ROOT::RVec<float> &t70_y,
    const ROOT::RVec<float> &q50_y,
    const ROOT::RVec<float> &q90_y)
{
    const double mref_mev_per_cm = 2.10;
    const double tail_weight_w = 1.0;

    const auto n = t70_u.size();
    ROOT::RVec<float> out(n, std::numeric_limits<float>::quiet_NaN());

    auto at = [](const ROOT::RVec<float> &v, std::size_t i) {
        return (i < v.size()) ? v[i] : std::numeric_limits<float>::quiet_NaN();
    };

    auto discriminant = [mref_mev_per_cm, tail_weight_w](float t70, float q50, float q90) {
        if (!(std::isfinite(t70) && std::isfinite(q50) && std::isfinite(q90)))
            return std::numeric_limits<double>::quiet_NaN();
        if (!(t70 > 0.0f && q50 > 0.0f && q90 > 0.0f))
            return std::numeric_limits<double>::quiet_NaN();
        return std::log(static_cast<double>(t70) / mref_mev_per_cm) +
               tail_weight_w * std::log(static_cast<double>(q90) / static_cast<double>(q50));
    };

    for (std::size_t i = 0; i < n; ++i)
    {
        std::vector<double> d_values;
        d_values.reserve(3);

        const double d_u = discriminant(at(t70_u, i), at(q50_u, i), at(q90_u, i));
        const double d_v = discriminant(at(t70_v, i), at(q50_v, i), at(q90_v, i));
        const double d_y = discriminant(at(t70_y, i), at(q50_y, i), at(q90_y, i));

        if (std::isfinite(d_u))
            d_values.push_back(d_u);
        if (std::isfinite(d_v))
            d_values.push_back(d_v);
        if (std::isfinite(d_y))
            d_values.push_back(d_y);

        if (d_values.empty())
            continue;

        std::sort(d_values.begin(), d_values.end());
        out[i] = static_cast<float>(std::exp(-d_values[d_values.size() / 2]));
    }

    return out;
}

ROOT::RVec<float> derive_track_mipness_logistic_median_plane_score(const ROOT::RVec<float> &mipness)
{
    const double d0 = 3.0;
    const double softness = 1.0;

    ROOT::RVec<float> out(mipness.size(), std::numeric_limits<float>::quiet_NaN());
    for (std::size_t i = 0; i < mipness.size(); ++i)
    {
        const double m = static_cast<double>(mipness[i]);
        if (!(std::isfinite(m) && m > 0.0))
            continue;
        const double d = -std::log(m);
        const double z = std::clamp((d - d0) / softness, -60.0, 60.0);
        out[i] = static_cast<float>(1.0 / (1.0 + std::exp(z)));
    }
    return out;
}

} // namespace


//____________________________________________________________________________
ROOT::RDF::RNode ColumnDerivationService::define(ROOT::RDF::RNode node, const ProcessorEntry &rec) const
{
    const bool is_data = (rec.source == Type::kData);
    const bool is_ext = (rec.source == Type::kExt);
    const bool is_mc = (rec.source == Type::kMC);

    const double scale_mc =
        (is_mc && rec.pot_nom > 0.0 && rec.pot_eqv > 0.0) ? (rec.pot_nom / rec.pot_eqv) : 1.0;
    const double scale_ext =
        (is_ext && rec.trig_nom > 0.0 && rec.trig_eqv > 0.0) ? (rec.trig_nom / rec.trig_eqv) : 1.0;

    node = node.Define("w_base", [is_mc, is_ext, scale_mc, scale_ext]() -> double {
        const double scale = is_mc ? scale_mc : (is_ext ? scale_ext : 1.0);
        return scale;
    });

    {
        const auto cnames = node.GetColumnNames();
        auto has = [&](const std::string &name) {
            return std::find(cnames.begin(), cnames.end(), name) != cnames.end();
        };

        if (!has("ppfx_cv"))
        {
            node = node.Define("ppfx_cv", [] { return 1.0f; });
        }
        if (!has("weightSpline"))
        {
            node = node.Define("weightSpline", [] { return 1.0f; });
        }
        if (!has("weightTune"))
        {
            node = node.Define("weightTune", [] { return 1.0f; });
        }
        if (!has("RootinoFix"))
        {
            node = node.Define("RootinoFix", [] { return 1.0; });
        }
    }

    {
        const auto cnames = node.GetColumnNames();
        auto has = [&](const std::string &name) {
            return std::find(cnames.begin(), cnames.end(), name) != cnames.end();
        };

        const bool has_mipness_median = has("track_mipness_median_plane_score");
        const bool has_mipness_inputs =
            has("track_dedx_T70_u") && has("track_dedx_Q50_u") && has("track_dedx_Q90_u") &&
            has("track_dedx_T70_v") && has("track_dedx_Q50_v") && has("track_dedx_Q90_v") &&
            has("track_dedx_T70_y") && has("track_dedx_Q50_y") && has("track_dedx_Q90_y");

        if (!has_mipness_median && has_mipness_inputs)
        {
            node = node.Define(
                "track_mipness_median_plane_score",
                derive_track_mipness_median_plane_score,
                {"track_dedx_T70_u", "track_dedx_Q50_u", "track_dedx_Q90_u",
                 "track_dedx_T70_v", "track_dedx_Q50_v", "track_dedx_Q90_v",
                 "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});
        }

        if (!has("track_mipness_logistic_median_plane_score") &&
            (has_mipness_median || has_mipness_inputs))
        {
            node = node.Define(
                "track_mipness_logistic_median_plane_score",
                derive_track_mipness_logistic_median_plane_score,
                {"track_mipness_median_plane_score"});
        }
    }

    if (is_mc)
    {
        node = node.Define(
            "w_nominal",
            [](double w_base, float w_spline, float w_tune, float w_flux_cv, double w_root) -> double {
                auto sanitise_weight = [](double w) {
                    if (!std::isfinite(w) || w <= 0.0)
                        return 1.0;
                    return w;
                };
                const double out = w_base *
                                   sanitise_weight(w_spline) *
                                   sanitise_weight(w_tune) *
                                   sanitise_weight(w_flux_cv) *
                                   sanitise_weight(w_root);
                if (!std::isfinite(out))
                    return 0.0;
                if (out < 0.0)
                    return 0.0;
                return out;
            },
            {"w_base", "weightSpline", "weightTune", "ppfx_cv", "RootinoFix"});
    }
    else
    {
        node = node.Define("w_nominal", [](double w) -> double { return w; }, {"w_base"});
    }


    if (is_mc)
    {
        node = node.Define(
            "in_fiducial",
            [](float x, float y, float z) {
                return SelectionService::is_in_truth_volume(x, y, z);
            },
            {"nu_vtx_x", "nu_vtx_y", "nu_vtx_z"});

        node = node.Define(
            "count_strange",
            [](int kplus, int kminus, int kzero, int lambda0, int sigplus, int sigzero, int sigminus) {
                return kplus + kminus + kzero + lambda0 + sigplus + sigzero + sigminus;
            },
            {"n_K_plus", "n_K_minus", "n_K0", "n_lambda", "n_sigma_plus", "n_sigma0", "n_sigma_minus"});

        node = node.Define(
            "is_strange",
            [](int strange) { return strange > 0; },
            {"count_strange"});

        node = node.Define(
            "analysis_channels",
            [](bool fv, int nu, int ccnc, int s, int np, int npim, int npip, int npi0, int ngamma) {
                const int npi = npim + npip;
                if (!fv)
                {
                    if (nu == 0)
                        return static_cast<int>(Channel::External);
                    return static_cast<int>(Channel::OutFV);
                }
                if (ccnc == 1)
                    return static_cast<int>(Channel::NC);
                if (ccnc == 0 && s > 0)
                {
                    if (s == 1)
                        return static_cast<int>(Channel::CCS1);
                    return static_cast<int>(Channel::CCSgt1);
                }
                if (std::abs(nu) == 12 && ccnc == 0)
                    return static_cast<int>(Channel::ECCC);
                if (std::abs(nu) == 14 && ccnc == 0)
                {
                    if (npi == 0 && np > 0)
                        return static_cast<int>(Channel::MuCC0pi_ge1p);
                    if (npi == 1 && npi0 == 0)
                        return static_cast<int>(Channel::MuCC1pi);
                    if (npi0 > 0 || ngamma >= 2)
                        return static_cast<int>(Channel::MuCCPi0OrGamma);
                    if (npi > 1)
                        return static_cast<int>(Channel::MuCCNpi);
                    return static_cast<int>(Channel::MuCCOther);
                }
                return static_cast<int>(Channel::Unknown);
            },
            {"in_fiducial", "nu_pdg", "int_ccnc", "count_strange", "n_p", "n_pi_minus", "n_pi_plus", "n_pi0", "n_gamma"});

        node = node.Define(
            "is_signal",
            [](bool is_nu_mu_cc, int ccnc, bool in_fiducial, float mu_p, float p_p, float pi_p, float lam_decay_sep) {
                const float min_mu_p = 0.10f;
                const float min_p_p = 0.30f;
                const float min_pi_p = 0.10f;
                const float min_lam_decay_sep = 0.50f;

                if (!is_nu_mu_cc)
                    return false;
                if (ccnc != 0)
                    return false;
                if (!in_fiducial)
                    return false;
                if (!std::isfinite(mu_p) || !std::isfinite(p_p) || !std::isfinite(pi_p) ||
                    !std::isfinite(lam_decay_sep))
                    return false;
                if (mu_p < min_mu_p || p_p < min_p_p || pi_p < min_pi_p)
                    return false;
                return lam_decay_sep >= min_lam_decay_sep;
            },
            {"is_nu_mu_cc", "int_ccnc", "in_fiducial", "mu_p", "p_p", "pi_p", "lam_decay_sep"});
    }
    else
    {
        const int nonmc_channel = is_ext ? static_cast<int>(Channel::External)
                                : (is_data ? static_cast<int>(Channel::DataInclusive)
                                : static_cast<int>(Channel::Unknown));

        node = node.Define("nu_vtx_x", [] { return -9999.0f; });
        node = node.Define("nu_vtx_y", [] { return -9999.0f; });
        node = node.Define("nu_vtx_z", [] { return -9999.0f; });

        node = node.Define("in_fiducial", [] { return false; });
        node = node.Define("is_strange", [] { return false; });
        node = node.Define("analysis_channels", [nonmc_channel] { return nonmc_channel; });
        node = node.Define("is_signal", [] { return false; });
        node = node.Define("recognised_signal", [] { return false; });
    }

    node = node.Define(
        "in_reco_fiducial",
        [](float x, float y, float z) {
            return SelectionService::is_in_reco_volume(x, y, z);
        },
        {"reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y", "reco_neutrino_vertex_sce_z"});

    {
        SelectionEntry srec{rec.source, Frame{node}};
        node = SelectionService::decorate(node, srec);
    }

    return node;
}
//____________________________________________________________________________

//____________________________________________________________________________
const ColumnDerivationService &ColumnDerivationService::instance()
{
    static const ColumnDerivationService ep{};
    return ep;
}
//____________________________________________________________________________
