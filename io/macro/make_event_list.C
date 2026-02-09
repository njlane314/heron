// io/macro/make_event_list.C
//
// Create a compact event-list ROOT file (merged) for fast plotting.
//
// Usage examples:
//   ./nuxsec macro make_event_list.C
//   ./nuxsec macro make_event_list.C 'make_event_list("scratch/out/event_list.root")'
//   ./nuxsec macro make_event_list.C 'make_event_list("scratch/out/event_list.root","/path/to/samples.tsv")'
//   ./nuxsec macro make_event_list.C 'make_event_list("scratch/out/event_list.root","/path/to/samples.tsv","true","reco_neutrino_vertex_sce_z,reco_neutrino_vertex_sce_x,reco_neutrino_vertex_sce_y")'
//   ./nuxsec macro make_event_list.C 'make_event_list("scratch/out/event_list.root","/path/to/samples.tsv","sel_muon","reco_neutrino_vertex_sce_z")'
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
//
// After this, you can plot from the event list using:
//   ROOT::RDataFrame("events", "scratch/out/event_list.root")
//
// Or modify stack_samples.C to detect ".root" input and use the event list directly.

#include <string>
#include <vector>

#include "EventListIO.hh"

static std::vector<std::string> default_event_columns()
{
    return {
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

        "reco_neutrino_vertex_sce_x",
        "reco_neutrino_vertex_sce_y",
        "reco_neutrino_vertex_sce_z",

        "in_reco_fiducial",
        "is_signal"
    };
}

int make_event_list(const std::string &out_root = "scratch/out/event_list.root",
                    const std::string &samples_tsv = "",
                    const std::string &base_sel = "true")
{
    return nu::EventListIO::build_event_list(out_root,
                                             samples_tsv,
                                             base_sel,
                                             default_event_columns());
}
