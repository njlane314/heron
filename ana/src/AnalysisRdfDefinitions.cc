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

constexpr double AnalysisRdfDefinitions::kRecognisedPurityMin = 0.5;
constexpr double AnalysisRdfDefinitions::kRecognisedCompletenessMin = 0.1;

constexpr float AnalysisRdfDefinitions::kTrainingFraction = 0.10f;
constexpr bool AnalysisRdfDefinitions::kTrainingIncludeExt = true;

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

    if (is_mc)
    {
        node = node.Define(
            "w_nominal",
            [](float w, float w_spline, float w_tune) {
                const float out = w * w_spline * w_tune;
                if (!std::isfinite(out))
                    return 0.0f;
                if (out < 0.0f)
                    return 0.0f;
                return out;
            },
            {"w_base", "weightSpline", "weightTune"});
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
            {"simb_mode"});

        node = node.Define(
            "is_cc",
            [](int interaction) {
                return interaction == 0;
            },
            {"simb_interaction"});

        node = node.Define(
            "is_nc",
            [](int interaction) {
                return interaction != 0;
            },
            {"simb_interaction"});

        node = node.Define(
            "is_ccnu",
            [](int interaction, int parent) {
                return interaction == 0 && parent == 0;
            },
            {"simb_interaction", "simb_mother"});

        node = node.Define(
            "is_ccnubar",
            [](int interaction, int parent) {
                return interaction == 0 && parent != 0;
            },
            {"simb_interaction", "simb_mother"});

        node = node.Define(
            "analysis_channels",
            [](int interaction,
               int mode,
               int parent,
               float x,
               float y,
               float z,
               int n_pi0,
               int n_piplus,
               int n_piminus,
               int n_photon,
               int n_proton,
               int n_muon,
               int n_electron,
               int n_neutron,
               int n_kaon,
               int n_other,
               float purity,
               float completeness,
               float vtx_x,
               float vtx_y,
               float vtx_z,
               int vtx_type,
               float slice_score,
               float cluster_fraction) {
                if (interaction != 0)
                {
                    return Channel::NC;
                }

                if (!is_in_truth_volume(x, y, z))
                {
                    return Channel::OutFV;
                }

                if (purity < kRecognisedPurityMin || completeness < kRecognisedCompletenessMin)
                {
                    return Channel::External;
                }

                const int n_pi = n_piplus + n_piminus + n_pi0;
                if (n_muon > 0)
                {
                    if (n_pi == 0)
                    {
                        return Channel::MuCC0pi_ge1p;
                    }
                    if (n_pi == 1)
                    {
                        return Channel::MuCC1pi;
                    }
                    if (n_pi > 1)
                    {
                        return Channel::MuCCNpi;
                    }
                    if (n_pi0 > 0 || n_photon > 0)
                    {
                        return Channel::MuCCPi0OrGamma;
                    }
                    return Channel::MuCCOther;
                }

                if (n_electron > 0)
                {
                    return Channel::ECCC;
                }

                if (n_pi == 0)
                {
                    return Channel::CCS1;
                }
                if (n_pi > 1)
                {
                    return Channel::CCSgt1;
                }

                return Channel::DataInclusive;
            },
            {"simb_interaction",
             "simb_mode",
             "simb_mother",
             "nu_vtx_x",
             "nu_vtx_y",
             "nu_vtx_z",
             "n_pi0",
             "n_piplus",
             "n_piminus",
             "n_photon",
             "n_proton",
             "n_muon",
             "n_electron",
             "n_neutron",
             "n_kaon",
             "n_other",
             "slice_purity",
             "slice_completeness",
             "reco_vtx_x",
             "reco_vtx_y",
             "reco_vtx_z",
             "reco_vtx_type",
             "topological_score",
             "slice_cluster_fraction"});
    }
    else
    {
        node = node.Define(
            "analysis_channels",
            []() {
                return Channel::DataInclusive;
            });
    }

    return node;
}

const AnalysisRdfDefinitions &AnalysisRdfDefinitions::Instance()
{
    static const AnalysisRdfDefinitions ep{};
    return ep;
}

}
