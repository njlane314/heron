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
