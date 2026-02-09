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

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>

#include "AnalysisConfigService.hh"
#include "ColumnDerivationService.hh"
#include "EventIO.hh"
#include "EventSampleFilterService.hh"
#include "RDataFrameService.hh"
#include "SampleCLI.hh"
#include "SampleIO.hh"

using namespace nu;

// ---- helpers ---------------------------------------------------------------

static inline std::string trim_copy(std::string s)
{
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

static std::vector<std::string> split_csv(const std::string &csv)
{
    std::vector<std::string> out;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        item = trim_copy(item);
        if (!item.empty())
            out.push_back(item);
    }
    return out;
}

static std::vector<std::string> default_event_columns()
{
    // Keep this compact and stable across data/ext/mc.
    // You can add more via extra_columns_csv argument.
    return {
        // Needed for stacking
        "analysis_channels",
        "w_nominal",

        // Useful for plotting selections without recomputing
        "sel_trigger",
        "sel_slice",
        "sel_fiducial",
        "sel_topology",
        "sel_muon",

        // Common extras (harmless and often useful)
        "sel_inclusive_mu_cc",
        "sel_reco_fv",
        "sel_triggered_slice",
        "sel_triggered_muon",

        // A few common reco vars (safe across data/ext/mc)
        "reco_neutrino_vertex_sce_x",
        "reco_neutrino_vertex_sce_y",
        "reco_neutrino_vertex_sce_z",

        // Optional but handy
        "in_reco_fiducial",
        "is_signal"
    };
}

static std::vector<std::string> merge_unique(std::vector<std::string> a, const std::vector<std::string> &b)
{
    a.insert(a.end(), b.begin(), b.end());
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    return a;
}

static void require_columns(const ROOT::RDF::RNode &node,
                            const std::vector<std::string> &cols,
                            const std::string &sample_name)
{
    const auto names = node.GetColumnNames();
    std::unordered_set<std::string> have(names.begin(), names.end());

    std::vector<std::string> missing;
    missing.reserve(cols.size());
    for (const auto &c : cols)
    {
        if (c == "sample_id")
            continue; // SnapshotService adds/defines this
        if (have.find(c) == have.end())
            missing.push_back(c);
    }

    if (!missing.empty())
    {
        std::ostringstream err;
        err << "make_event_list: missing columns after derivation for sample '" << sample_name << "':\n";
        for (const auto &m : missing)
            err << "  - " << m << "\n";
        err << "Fix: remove these from extra_columns_csv, or ensure they exist/are defined for all sample types.\n";
        throw std::runtime_error(err.str());
    }
}

// ---- main macro ------------------------------------------------------------

int make_event_list(const std::string &out_root = "scratch/out/event_list.root",
                    const std::string &samples_tsv = "",
                    const std::string &base_sel = "true",
                    const std::string &extra_columns_csv = "",
                    bool apply_origin_filters = true,
                    bool enable_mt = true)
{
    if (enable_mt)
        ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_samples_tsv() : samples_tsv;
    std::cout << "[make_event_list] samples_tsv=" << list_path << "\n";
    std::cout << "[make_event_list] out_root=" << out_root << "\n";
    std::cout << "[make_event_list] base_sel=" << base_sel << "\n";

    const auto sample_list = read_samples(list_path);

    const auto &analysis = AnalysisConfigService::instance();
    const std::string tree_name = analysis.tree_name();

    // Columns to snapshot
    std::vector<std::string> cols = default_event_columns();
    if (!extra_columns_csv.empty())
    {
        const auto extra = split_csv(extra_columns_csv);
        cols = merge_unique(std::move(cols), extra);
    }

    // Build sample_refs metadata (for sample_id -> sample bookkeeping)
    std::vector<SampleInfo> refs;
    refs.reserve(sample_list.size());

    for (const auto &sl : sample_list)
    {
        SampleIO::Sample s = SampleIO::read(sl.output_path);

        SampleInfo si;
        si.sample_name = s.sample_name;
        si.sample_rootio_path = sl.output_path;
        si.sample_origin = static_cast<int>(s.origin);
        si.beam_mode = static_cast<int>(s.beam);
        si.subrun_pot_sum = s.subrun_pot_sum;
        si.db_tortgt_pot_sum = s.db_tortgt_pot_sum;
        si.db_tor101_pot_sum = s.db_tor101_pot_sum;

        refs.push_back(std::move(si));
    }

    Header header;
    header.analysis_name = analysis.name();
    header.provenance_tree = tree_name;   // original analysis tree
    header.event_tree = "events";         // merged event list tree
    header.sample_list_source = list_path;

    // Store the schema (as a plain string) into the file for provenance/debug
    std::ostringstream schema;
    schema << "# nuxsec event list columns (macro make_event_list.C)\n";
    for (const auto &c : cols)
        schema << c << "\n";

    // Overwrite output file and write header + sample_refs
    EventIO::init(out_root, header, refs, schema.str(), "plot");

    // Append merged event tree
    EventIO out(out_root, EventIO::OpenMode::kUpdate);

    for (size_t i = 0; i < sample_list.size(); ++i)
    {
        SampleIO::Sample sample = SampleIO::read(sample_list[i].output_path);

        std::cout << "[make_event_list] sample=" << sample.sample_name
                  << " id=" << i
                  << " origin=" << SampleIO::sample_origin_name(sample.origin)
                  << " beam=" << SampleIO::beam_mode_name(sample.beam)
                  << "\n";

        ROOT::RDataFrame rdf = RDataFrameService::load_sample(sample, tree_name);

        const ProcessorEntry proc = analysis.make_processor(sample);

        ROOT::RDF::RNode node = ColumnDerivationService::instance().define(rdf, proc);

        if (apply_origin_filters)
        {
            node = EventSampleFilterService::apply(std::move(node), sample.origin);
        }

        // Ensure the event list schema is consistent and exists for this sample
        require_columns(node, cols, sample.sample_name);

        // Snapshot into ONE merged tree called "events"
        out.snapshot_event_list_merged(std::move(node),
                                       static_cast<int>(i),
                                       sample.sample_name,
                                       cols,
                                       base_sel,
                                       header.event_tree);
    }

    std::cout << "[make_event_list] done: " << out_root << "\n";
    return 0;
}
