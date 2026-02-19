/* -- C++ -- */
/**
 *  @file  ana/src/ColumnDerivationService.cpp
 *
 *  @brief Variable definitions for analysis RDataFrame processing.
 */

#include "ColumnDerivationService.hh"

#include <algorithm>
#include <cmath>
#include <string>

#include <ROOT/RVec.hxx>

#include "SelectionService.hh"

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

        {
            const auto mc_cnames = node.GetColumnNames();
            auto has_mc = [&](const std::string &name) {
                return std::find(mc_cnames.begin(), mc_cnames.end(), name) != mc_cnames.end();
            };
            if (!has_mc("interaction_mode"))
            {
                if (has_mc("int_mode"))
                {
                    node = node.Define("interaction_mode", [](int m) { return m; }, {"int_mode"});
                }
                else
                {
                    node = node.Define("interaction_mode", [] { return -1; });
                }
            }
        }

        node = node.Define(
            "analysis_channels",
            [](bool in_fiducial,
               int nu_pdg,
               int ccnc,
               int interaction_mode,
               int n_p,
               int n_pi_minus,
               int n_pi_plus,
               int n_pi0,
               int n_gamma,
               bool is_nu_mu_cc,
               float mu_p,
               float p_p,
               float pi_p,
               float lam_decay_sep) {
                return AnalysisChannels::to_int(
                    AnalysisChannels::classify_analysis_channel(
                        in_fiducial,
                        nu_pdg,
                        ccnc,
                        interaction_mode,
                        n_p,
                        n_pi_minus,
                        n_pi_plus,
                        n_pi0,
                        n_gamma,
                        is_nu_mu_cc,
                        mu_p,
                        p_p,
                        pi_p,
                        lam_decay_sep));
            },
            {"in_fiducial",
             "nu_pdg",
             "int_ccnc",
             "interaction_mode",
             "n_p",
             "n_pi_minus",
             "n_pi_plus",
             "n_pi0",
             "n_gamma",
             "is_nu_mu_cc",
             "mu_p",
             "p_p",
             "pi_p",
             "lam_decay_sep"});


        node = node.Define(
            "is_signal",
            [](bool is_nu_mu_cc, int ccnc, bool in_fiducial, float mu_p, float p_p, float pi_p, float lam_decay_sep) {
                return AnalysisChannels::is_signal(
                    is_nu_mu_cc,
                    ccnc,
                    in_fiducial,
                    mu_p,
                    p_p,
                    pi_p,
                    lam_decay_sep);
            },
            {"is_nu_mu_cc", "int_ccnc", "in_fiducial", "mu_p", "p_p", "pi_p", "lam_decay_sep"});
    }
    else
    {
        const int nonmc_channel = is_ext ? static_cast<int>(AnalysisChannels::AnalysisChannel::External)
                                : (is_data ? static_cast<int>(AnalysisChannels::AnalysisChannel::DataInclusive)
                                : static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown));

        node = node.Define("nu_vtx_x", [] { return -9999.0f; });
        node = node.Define("nu_vtx_y", [] { return -9999.0f; });
        node = node.Define("nu_vtx_z", [] { return -9999.0f; });

        node = node.Define("in_fiducial", [] { return false; });
        node = node.Define("is_strange", [] { return false; });
        node = node.Define("analysis_channels", [nonmc_channel] { return nonmc_channel; });
        node = node.Define("interaction_mode", [] { return -1; });
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
