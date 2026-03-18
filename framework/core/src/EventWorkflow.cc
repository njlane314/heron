/* -- C++ -- */
/**
 *  @file  framework/core/src/EventWorkflow.cc
 *
 *  @brief Event-level output builder (invoked by the unified heron CLI).
 */

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RVec.hxx>

#include "AnalysisConfigService.hh"
#include "AppUtils.hh"
#include "ColumnDerivationService.hh"
#include "Dataset.hh"
#include "EventCLI.hh"
#include "EventColumnProvider.hh"
#include "EventListIO.hh"
#include "EventSampleFilterService.hh"
#include "RDataFrameService.hh"
#include "StatusMonitor.hh"

namespace
{

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto names = node.GetColumnNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

template <class T>
ROOT::RDF::RNode define_if_missing(ROOT::RDF::RNode node,
                                   const std::string &name,
                                   const T &value)
{
    if (has_column(node, name))
    {
        return node;
    }

    return node.Define(name, [value]() { return value; });
}

ROOT::RDF::RNode add_event_weight_defaults(ROOT::RDF::RNode node)
{
    using UShortVec = ROOT::VecOps::RVec<unsigned short>;
    constexpr unsigned short kPackedUnity = 1000;

    // Neutral defaults for samples with no event-weight branches (e.g. EXT/data).
    node = define_if_missing(node, "weightSpline", 1.f);
    node = define_if_missing(node, "weightTune", 1.f);
    node = define_if_missing(node, "weightSplineTimesTune", 1.f);
    node = define_if_missing(node, "ppfx_cv", 1.f);

    for (const char *name : {
             "knobRPAup", "knobRPAdn",
             "knobCCMECup", "knobCCMECdn",
             "knobAxFFCCQEup", "knobAxFFCCQEdn",
             "knobVecFFCCQEup", "knobVecFFCCQEdn",
             "knobDecayAngMECup", "knobDecayAngMECdn",
             "knobThetaDelta2Npiup", "knobThetaDelta2Npidn",
             "knobThetaDelta2NRadup", "knobThetaDelta2NRaddn",
             "knobNormCCCOHup", "knobNormCCCOHdn",
             "knobNormNCCOHup", "knobNormNCCOHdn",
             "knobxsr_scc_Fv3up", "knobxsr_scc_Fv3dn",
             "knobxsr_scc_Fa3up", "knobxsr_scc_Fa3dn",
             "RootinoFix"})
    {
        node = define_if_missing(node, name, 1.0);
    }

    // Packed vectors store weight * 1000 in unsigned short.
    node = define_if_missing(node, "weightsGenie", UShortVec(500, kPackedUnity));
    node = define_if_missing(node, "weightsPPFX", UShortVec(600, kPackedUnity));

    // These sizes are production-dependent, so empty is the safest default.
    node = define_if_missing(node, "weightsFlux", UShortVec{});
    node = define_if_missing(node, "weightsReint", UShortVec{});
    node = define_if_missing(node, "weightsGenieUp", UShortVec{});
    node = define_if_missing(node, "weightsGenieDn", UShortVec{});

    return node;
}

} // namespace

int run(const EventArgs &event_args, const std::string &log_prefix)
{
    ROOT::EnableImplicitMT();

    const auto &analysis = AnalysisConfigService::instance();
    const Dataset dataset = Dataset::load(event_args.list_path);

    const auto start_time = std::chrono::steady_clock::now();
    log_event_start(log_prefix, dataset.inputs().size());

    StatusMonitor status_monitor(
        log_prefix,
        "action=event_build status=running message=processing");

    const auto &inputs = dataset.inputs();
    const auto &sample_infos = dataset.sample_infos();

    if (event_args.columns_tsv_path.empty())
    {
        throw std::runtime_error(
            "Event columns TSV is required; pass selection (use 'true' for empty selection) and COLUMNS.tsv");
    }

    const std::string columns_tsv_path = event_args.columns_tsv_path;

    const std::string provenance_tree = "heron_art_provenance/run_subrun";
    const std::string event_tree = analysis.tree_name();
    const std::string output_event_tree = "events";

    nu::EventListHeader header;
    header.analysis_name = analysis.name();
    header.provenance_tree = provenance_tree;
    header.event_tree = output_event_tree;
    header.sample_list_source = event_args.list_path;
    header.heron_set = workspace_set();

    const std::filesystem::path output_path(event_args.output_root);
    const auto output_parent = output_path.parent_path();
    if (!output_parent.empty())
    {
        header.event_output_dir = output_parent.string();
        std::filesystem::create_directories(output_parent);
    }

    const EventColumnProvider column_provider(columns_tsv_path);

    nu::EventListIO::init(event_args.output_root,
                          header,
                          sample_infos,
                          column_provider.schema_tsv(),
                          column_provider.schema_tag());
    nu::EventListIO event_io(event_args.output_root,
                             nu::EventListIO::OpenMode::kUpdate);

    for (size_t i = 0; i < inputs.size(); ++i)
    {
        const auto &input = inputs[i];
        const SampleIO::Sample &sample = input.sample;
        const int sample_id = static_cast<int>(i);

        log_stage(
            log_prefix,
            "ensure_tree",
            "sample=" + sample.sample_name + " tree=" + event_tree);

        ensure_tree_present(sample, event_tree);

        log_stage(
            log_prefix,
            "load_rdf",
            "sample=" + sample.sample_name);

        ROOT::RDataFrame rdf = RDataFrameService::load_sample(sample, event_tree);

        log_stage(
            log_prefix,
            "make_processor",
            "sample=" + sample.sample_name);

        const ProcessorEntry proc_entry = analysis.make_processor(sample);
        const auto &processor = ColumnDerivationService::instance();

        log_stage(
            log_prefix,
            "define_columns",
            "sample=" + sample.sample_name);

        ROOT::RDF::RNode node = processor.define(rdf, proc_entry);
        node = add_event_weight_defaults(node);

        const char *filter_stage = EventSampleFilterService::filter_stage(sample.origin);
        if (filter_stage != nullptr)
        {
            log_stage(
                log_prefix,
                filter_stage,
                "sample=" + sample.sample_name);
            node = EventSampleFilterService::apply(node, sample.origin);
        }

        std::string snapshot_message = "sample=" + sample.sample_name;
        if (!event_args.selection.empty())
        {
            snapshot_message += " selection=" + event_args.selection;
        }
        log_stage(
            log_prefix,
            "snapshot",
            snapshot_message);

        const ULong64_t n_written =
            event_io.snapshot_event_list_merged(node,
                                                sample_id,
                                                sample.sample_name,
                                                column_provider.columns(),
                                                event_args.selection,
                                                output_event_tree);

        std::ostringstream log_message;
        log_message << "action=event_snapshot status=complete analysis=" << analysis.name()
                    << " sample=" << sample.sample_name
                    << " kind=" << SampleIO::sample_origin_name(sample.origin)
                    << " beam=" << SampleIO::beam_mode_name(sample.beam)
                    << " events_written=" << n_written
                    << " output=" << event_args.output_root;
        if (!event_args.selection.empty())
        {
            log_message << " selection=" << event_args.selection;
        }
        log_success(log_prefix, log_message.str());
    }
    status_monitor.stop();

    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    log_event_finish(log_prefix, inputs.size(), elapsed_seconds);

    return 0;
}
