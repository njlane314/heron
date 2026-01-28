/* -- C++ -- */
/**
 *  @file  ana/src/ColumnDerivationService.cc
 *
 *  @brief Variable definitions for analysis RDataFrame processing.
 */

#include "ColumnDerivationService.hh"

#include <algorithm>
#include <cmath>
#include <string>

#include <ROOT/RVec.hxx>

#include "Selection.hh"

namespace nuxsec
{

//____________________________________________________________________________
ROOT::RDF::RNode ColumnDerivationService::define(ROOT::RDF::RNode node, const ProcessorEntry &rec) const
{
    const bool is_data = (rec.source == SourceKind::kData);
    const bool is_ext = (rec.source == SourceKind::kExt);
    const bool is_mc = (rec.source == SourceKind::kMC);

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
                return selection::SelectionService::is_in_truth_volume(x, y, z);
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
            {"in_fiducial", "nu_pdg", "int_ccnc", "count_strange", "n_p", "n_pi_minus", "n_pi_plus", "n_pi0", "n_gamma"});

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
    }
    else
    {
        const int nonmc_channel = is_ext ? static_cast<int>(Channel::External)
                                : (is_data ? static_cast<int>(Channel::DataInclusive)
                                : static_cast<int>(Channel::Unknown));

        node = node.Define("in_fiducial", [] { return false; });
        node = node.Define("is_strange", [] { return false; });
        node = node.Define("analysis_channels", [nonmc_channel] { return nonmc_channel; });
        node = node.Define("is_signal", [] { return false; });
        node = node.Define("recognised_signal", [] { return false; });
    }

    node = node.Define(
        "in_reco_fiducial",
        [](float x, float y, float z) {
            return selection::SelectionService::is_in_reco_volume(x, y, z);
        },
        {"reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y", "reco_neutrino_vertex_sce_z"});

    node = node.Define(
        "sel_trigger",
        [is_mc](float pe_beam, float pe_veto, int sw) {
            const bool requires_dataset_gate = is_mc;
            const bool dataset_gate = requires_dataset_gate
                                          ? (pe_beam > selection::SelectionService::trigger_min_beam_pe &&
                                             pe_veto < selection::SelectionService::trigger_max_veto_pe && sw > 0)
                                          : true;
            return dataset_gate;
        },
        {"optical_filter_pe_beam", "optical_filter_pe_veto", "software_trigger"});

    node = node.Define(
        "sel_slice",
        [](bool trigger, int ns, float topo) {
            if (!trigger)
                return false;
            return ns == selection::SelectionService::slice_required_count &&
                   topo > selection::SelectionService::slice_min_topology_score;
        },
        {"sel_trigger", "num_slices", "topological_score"});

    node = node.Define(
        "sel_fiducial",
        [](bool slice, bool fv) {
            if (!slice)
                return false;
            return fv;
        },
        {"sel_slice", "in_reco_fiducial"});

    node = node.Define(
        "sel_topology",
        [](bool fiducial, float cf, float cl) {
            if (!fiducial)
                return false;
            return cf >= selection::SelectionService::topology_min_contained_fraction &&
                   cl >= selection::SelectionService::topology_min_cluster_fraction;
        },
        {"sel_fiducial", "contained_fraction", "slice_cluster_fraction"});

    node = node.Define(
        "sel_muon",
        [](bool topology,
           const ROOT::RVec<float> &scores,
           const ROOT::RVec<float> &lengths,
           const ROOT::RVec<float> &distances,
           const ROOT::RVec<unsigned> &generations) {
            if (!topology)
                return false;
            const auto n = scores.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                const bool passes = scores[i] > selection::SelectionService::muon_min_track_score &&
                                    lengths[i] > selection::SelectionService::muon_min_track_length &&
                                    distances[i] < selection::SelectionService::muon_max_track_distance &&
                                    generations[i] == selection::SelectionService::muon_required_generation;
                if (passes)
                {
                    return true;
                }
            }
            return false;
        },
        {"sel_topology",
         "track_shower_scores",
         "track_length",
         "track_distance_to_vertex",
         "pfp_generations"});

    node = node.Define(
        "sel_inclusive_mu_cc",
        [](bool muon) { return muon; },
        {"sel_muon"});

    node = node.Define(
        "sel_reco_fv",
        [](bool reco_fv) { return reco_fv; },
        {"in_reco_fiducial"});

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

} // namespace nuxsec
