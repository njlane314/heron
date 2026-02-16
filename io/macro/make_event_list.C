// io/macro/make_event_list.C
//
// Create a compact event-list ROOT file (merged) for fast plotting.
//
// Usage examples:
//   ./heron macro make_event_list.C
//   ./heron macro make_event_list.C 'make_event_list("./scratch/out/event_list_myana.root")'
//   ./heron macro make_event_list.C 'make_event_list("./scratch/out/event_list_myana.root","/path/to/samples.tsv")'
//   ./heron macro make_event_list.C 'make_event_list("./scratch/out/event_list_myana.root","/path/to/samples.tsv","true","reco_neutrino_vertex_sce_z,reco_neutrino_vertex_sce_x,reco_neutrino_vertex_sce_y")'
//   ./heron macro make_event_list.C 'make_event_list("./scratch/out/event_list_myana.root","/path/to/samples.tsv","sel_muon","reco_neutrino_vertex_sce_z")'
//
// What it writes:
//   - TObjString keys: analysis_name, provenance_tree, event_tree, sample_list_source
//   - TTree "sample_refs": sample_id -> (sample_name, origin, beam, POT sums, etc.)
//   - TTree "events": merged event list (only requested columns + auto "sample_id")
//
// Notes:
//   - This runs ColumnDerivationService once per sample, then snapshots the derived columns.
//   - Output file is overwritten (RECREATE).
//   - base_sel is applied during snapshot to reduce the event list size if desired.
//   - extra_columns_csv lets you add plot variables beyond the defaults.
//   - Default output path is /exp/uboone/data/users/$USER/event_list_<analysis>.root.
//   - USER must be set when out_root is not provided.
//
// After this, you can plot from the event list using:
//   ROOT::RDataFrame("events", "./scratch/out/event_list_myana.root")
//
// Or modify stack_samples.C to detect ".root" input and use the event list directly.

#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AnalysisConfigService.hh"
#include "EventListIO.hh"

static std::vector<std::string> default_event_columns()
{
    return {
        "run",
        "sub",
        "evt",

        "analysis_channels",
        "w_nominal",

        "sel_trigger",
        "sel_slice",
        "sel_fiducial",
        "sel_topology",
        "sel_muon",

        "sel_inclusive_mu_cc",
        "sel_reco_fv",
        "sel_triggered_slice",
        "sel_triggered_muon",

        "optical_filter_pe_beam",
        "optical_filter_pe_veto",
        "software_trigger",
        "num_slices",
        "topological_score",
        "track_shower_scores",
        "track_length",
        "track_distance_to_vertex",
        "track_mipness_median_plane_score",
        "backtracked_pdg_codes",
        "pfp_generations",
        "track_dedx_T70_u",
        "track_dedx_T70_v",
        "track_dedx_T70_y",
        "track_dedx_Q50_u",
        "track_dedx_Q50_v",
        "track_dedx_Q50_y",
        "track_dedx_Q90_u",
        "track_dedx_Q90_v",
        "track_dedx_Q90_y",

        "reco_neutrino_vertex_sce_x",
        "reco_neutrino_vertex_sce_y",
        "reco_neutrino_vertex_sce_z",
        "is_vtx_in_image_u",
        "is_vtx_in_image_v",
        "is_vtx_in_image_w",

        "detector_image_u",
        "detector_image_v",
        "detector_image_w",
        "semantic_image_u",
        "semantic_image_v",
        "semantic_image_w",
        "active_pixels_u",
        "active_pixels_v",
        "active_pixels_w",
        "slice_semantic_active_pixels_u",
        "slice_semantic_active_pixels_v",
        "slice_semantic_active_pixels_w",
        "inf_scores",

        "nu_vtx_x",
        "nu_vtx_y",
        "nu_vtx_z",

        "in_reco_fiducial",
        "is_signal",

        // --------------------------------------------------------------------
        // Truth neutrino kinematics (TruthAnalysis)
        // --------------------------------------------------------------------
        "nu_pdg",
        "int_ccnc",
        "int_mode",
        "int_type",
        "is_nu_mu_cc",
        "nu_vtx_in_fv",

        "nu_E",
        "nu_theta",
        "nu_pt",

        "kin_W",
        "bjorken_x",
        "inelasticity_y",
        "Q2",

        // --------------------------------------------------------------------
        // Lambda truth + decay kinematics (LambdaAnalysis_tool)
        // --------------------------------------------------------------------
        "mu_truth_trackid",
        "mu_truth_pdg",
        "mu_p",
        "mu_theta",

        "lam_trackid",
        "lam_pdg",
        "lam_E",
        "lam_mass",
        "lam_p_mag",
        "lam_ct",
        "lam_decay_sep",

        "p_trackid",
        "pi_trackid",
        "p_p",
        "pi_p",
        "ppi_opening_angle",

        // --------------------------------------------------------------------
        // Pattern-recognition metrics (LambdaAnalysis_tool / slice-level)
        // --------------------------------------------------------------------
        "pr_valid_assignment",
        "pr_mu_purity",
        "pr_mu_completeness",
        "pr_p_purity",
        "pr_p_completeness",
        "pr_pi_purity",
        "pr_pi_completeness",

        // Truth-hit content per slice (useful for defining the efficiency denom)
        "mu_true_hits",
        "p_true_hits",
        "pi_true_hits"
    };
}

int make_event_list(const std::string &out_root = "",
                    const std::string &samples_tsv = "",
                    const std::string &base_sel = "true")
{
    std::string out_path = out_root;
    if (out_path.empty())
    {
        const auto &analysis = AnalysisConfigService::instance();
        std::ostringstream name;
        const char *user = std::getenv("USER");
        if (!user || !*user)
            throw std::runtime_error("make_event_list: USER must be set when out_root is empty");
        name << "/exp/uboone/data/users/" << user << "/event_list_" << analysis.name() << ".root";
        out_path = name.str();
    }

    return nu::EventListIO::build_event_list(out_path,
                                             samples_tsv,
                                             base_sel,
                                             default_event_columns());
}
