/* -- C++ -- */
/**
 *  @file  ana/src/AnalysisRdfDefinitions.cc
 *
 *  @brief Variable definitions for analysis RDataFrame processing.
 */

#include "AnalysisRdfDefinitions.hh"

#include <algorithm>
#include <cmath>
#include <string>

#include <ROOT/RVec.hxx>

namespace nuxsec
{

namespace
{

constexpr float min_x = 5.f;
constexpr float max_x = 251.f;
constexpr float min_y = -110.f;
constexpr float max_y = 110.f;
constexpr float min_z = 20.f;
constexpr float max_z = 986.f;

constexpr float reco_gap_min_z = 675.f;
constexpr float reco_gap_max_z = 775.f;

template <typename T>
inline bool is_within(const T &value, float low, float high)
{
    return value > low && value < high;
}

template <typename X, typename Y, typename Z>
inline bool is_in_active_volume(const X &x, const Y &y, const Z &z)
{
    return is_within(x, min_x, max_x) &&
           is_within(y, min_y, max_y) &&
           is_within(z, min_z, max_z);
}

template <typename X, typename Y, typename Z>
inline bool is_in_truth_volume(const X &x, const Y &y, const Z &z)
{
    return is_in_active_volume(x, y, z);
}

template <typename X, typename Y, typename Z>
inline bool is_in_reco_volume(const X &x, const Y &y, const Z &z)
{
    return is_in_active_volume(x, y, z) && (z < reco_gap_min_z || z > reco_gap_max_z);
}

}

const double AnalysisRdfDefinitions::kRecognisedPurityMin = 0.5;
const double AnalysisRdfDefinitions::kRecognisedCompletenessMin = 0.1;

const float AnalysisRdfDefinitions::kTrainingFraction = 0.10f;
const bool AnalysisRdfDefinitions::kTrainingIncludeExt = true;

bool AnalysisRdfDefinitions::IsInTruthVolume(float x, float y, float z) noexcept
{
    return is_in_truth_volume(x, y, z);
}

bool AnalysisRdfDefinitions::IsInRecoVolume(float x, float y, float z) noexcept
{
    return is_in_reco_volume(x, y, z);
}

//____________________________________________________________________________
ROOT::RDF::RNode AnalysisRdfDefinitions::Define(ROOT::RDF::RNode node, const ProcessorEntry &rec) const
{
    const bool is_data = (rec.source == SourceKind::kData);
    const bool is_ext = (rec.source == SourceKind::kExt);
    const bool is_mc = (rec.source == SourceKind::kMC);

    const double scale_mc =
        (is_mc && rec.pot_nom > 0.0 && rec.pot_eqv > 0.0) ? (rec.pot_nom / rec.pot_eqv) : 1.0;
    const double scale_ext =
        (is_ext && rec.trig_nom > 0.0 && rec.trig_eqv > 0.0) ? (rec.trig_nom / rec.trig_eqv) : 1.0;

    node = node.Define("w_base", [is_mc, is_ext, scale_mc, scale_ext] {
        const double scale = is_mc ? scale_mc : (is_ext ? scale_ext : 1.0);
        return static_cast<float>(scale);
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
        if (!has("RootinoFix"))
        {
            node = node.Define("RootinoFix", [] { return 1.0; });
        }
    }

    if (is_mc)
    {
        node = node.Define(
            "w_nominal",
            [](float w_base, float w_spline, float w_tune, float w_flux_cv, double w_root) {
                const double out = static_cast<double>(w_base) *
                                   static_cast<double>(w_spline) *
                                   static_cast<double>(w_tune) *
                                   static_cast<double>(w_flux_cv) *
                                   static_cast<double>(w_root);
                if (!std::isfinite(out))
                    return 0.0f;
                if (out < 0.0)
                    return 0.0f;
                return static_cast<float>(out);
            },
            {"w_base", "weightSpline", "weightTune", "ppfx_cv", "RootinoFix"});
    }
    else
    {
        node = node.Define("w_nominal", [](float w) { return w; }, {"w_base"});
    }

    {
        const bool trainable = is_mc || (is_ext && kTrainingIncludeExt);

        const auto cnames = node.GetColumnNames();
        auto has = [&](const std::string &name) {
            return std::find(cnames.begin(), cnames.end(), name) != cnames.end();
        };

        const bool have_ml_u = has("ml_u");

        if (!have_ml_u)
        {
            node = node.Define("ml_u", [] { return 0.0f; });
        }

        if (!has("is_training"))
        {
            node = node.Define(
                "is_training",
                [trainable, have_ml_u](float u) {
                    if (!trainable || !have_ml_u)
                        return false;
                    return u < kTrainingFraction;
                },
                {"ml_u"});
        }

        if (!has("is_template"))
        {
            node = node.Define(
                "is_template",
                [trainable](bool t) { return !trainable || !t; },
                {"is_training"});
        }

        if (!has("w_template"))
        {
            node = node.Define(
                "w_template",
                [trainable, have_ml_u](float w, bool t) {
                    if (!trainable || !have_ml_u)
                        return w;
                    if (t)
                        return 0.0f;
                    const float keep = 1.0f - kTrainingFraction;
                    if (keep <= 0.0f)
                        return 0.0f;
                    return w / keep;
                },
                {"w_nominal", "is_training"});
        }
    }

    if (is_mc)
    {
        node = node.Define(
            "in_fiducial",
            [](float x, float y, float z) {
                return IsInTruthVolume(x, y, z);
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
            "scattering_mode",
            [](int mode) {
                switch (mode)
                {
                case 0:
                    return 0;
                case 1:
                    return 1;
                case 2:
                    return 2;
                case 3:
                    return 3;
                case 10:
                    return 10;
                default:
                    return -1;
                }
            },
            {"int_mode"});

        node = node.Define(
            "analysis_channels",
            [](bool fv, int nu, int ccnc, int s, int np, int npim, int npip, int npi0, int ngamma) {
                const int npi = npim + npip;
                if (!fv)
                {
                    if (nu == 0)
                        return static_cast<int>(Channel::OutFV);
                    return static_cast<int>(Channel::External);
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
            {"in_fiducial", "nu_pdg", "int_ccnc", "count_strange",
             "n_p", "n_pi_minus", "n_pi_plus", "n_pi0", "n_gamma"});

        node = node.Define(
            "is_signal",
            [](bool is_nu_mu_cc, const ROOT::RVec<int> &lambda_decay_in_fid) {
                if (!is_nu_mu_cc)
                    return false;
                for (auto v : lambda_decay_in_fid)
                    if (v)
                        return true;
                return false;
            },
            {"is_nu_mu_cc", "lambda_decay_in_fid"});

        node = node.Define(
            "recognised_signal",
            [](bool is_sig, float purity, float completeness) {
                return is_sig && purity > static_cast<float>(kRecognisedPurityMin) &&
                       completeness > static_cast<float>(kRecognisedCompletenessMin);
            },
            {"is_signal", "neutrino_purity_from_pfp", "neutrino_completeness_from_pfp"});
    }
    else
    {
        const int nonmc_channel =
            is_ext ? static_cast<int>(Channel::External)
                   : (is_data ? static_cast<int>(Channel::DataInclusive)
                              : static_cast<int>(Channel::Unknown));

        node = node.Define("in_fiducial", [] { return false; });
        node = node.Define("is_strange", [] { return false; });
        node = node.Define("scattering_mode", [] { return -1; });
        node = node.Define("analysis_channels", [nonmc_channel] { return nonmc_channel; });
        node = node.Define("is_signal", [] { return false; });
        node = node.Define("recognised_signal", [] { return false; });
    }

    node = node.Define(
        "in_reco_fiducial",
        [](float x, float y, float z) {
            return IsInRecoVolume(x, y, z);
        },
        {"reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y", "reco_neutrino_vertex_sce_z"});

    node = node.Define(
        "sel_template",
        [](bool is_template) { return is_template; },
        {"is_template"});

    node = node.Define(
        "sel_reco_fv",
        [](bool is_template, bool reco_fv) { return is_template && reco_fv; },
        {"is_template", "in_reco_fiducial"});

    node = node.Define(
        "sel_signal",
        [](bool is_template, bool reco_fv, bool recognised_signal) {
            return is_template && reco_fv && recognised_signal;
        },
        {"is_template", "in_reco_fiducial", "recognised_signal"});

    node = node.Define(
        "sel_bkg",
        [](bool is_template, bool reco_fv, bool recognised_signal) {
            return is_template && reco_fv && !recognised_signal;
        },
        {"is_template", "in_reco_fiducial", "recognised_signal"});

    return node;
}
//____________________________________________________________________________

//____________________________________________________________________________
const AnalysisRdfDefinitions &AnalysisRdfDefinitions::Instance()
{
    static const AnalysisRdfDefinitions ep{};
    return ep;
}
//____________________________________________________________________________

} // namespace nuxsec
