/* -- C++ -- */
/**
 *  @file  framework/core/src/EventWorkflow.cc
 *
 *  @brief Event-level output builder (invoked by the unified heron CLI).
 */

#include <chrono>
#include <filesystem>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "AppUtils.hh"
#include "EventCLI.hh"
#include "EventColumnProvider.hh"
#include "EventSampleFilterService.hh"
#include "RDataFrameService.hh"
#include "StatusMonitor.hh"

int run(const EventArgs &event_args, const std::string &log_prefix)
{
    ROOT::EnableImplicitMT();

    const auto &analysis = AnalysisConfigService::instance();
    const auto entries = read_samples(event_args.list_path);

    const auto start_time = std::chrono::steady_clock::now();
    log_event_start(log_prefix, entries.size());

    StatusMonitor status_monitor(
        log_prefix,
        "action=event_build status=running message=processing");

    std::vector<EventInput> inputs;
    inputs.reserve(entries.size());

    std::vector<nu::SampleInfo> sample_infos;
    sample_infos.reserve(entries.size());

    for (const auto &entry : entries)
    {
        SampleIO::Sample sample = SampleIO::read(entry.output_path);

        nu::SampleInfo info;
        info.sample_name = sample.sample_name;
        info.sample_rootio_path = entry.output_path;
        info.sample_origin = static_cast<int>(sample.origin);
        info.beam_mode = static_cast<int>(sample.beam);
        info.subrun_pot_sum = sample.subrun_pot_sum;
        info.db_tortgt_pot_sum = sample.db_tortgt_pot_sum;
        info.db_tor101_pot_sum = sample.db_tor101_pot_sum;
        sample_infos.push_back(std::move(info));

        EventInput input;
        input.entry = entry;
        input.sample = std::move(sample);
        inputs.push_back(std::move(input));
    }

    const std::filesystem::path default_columns_tsv =
        repo_root_dir() / "core" / "config" / "event_columns.tsv";
    const std::string columns_tsv_path = event_args.columns_tsv_path.empty()
                                             ? default_columns_tsv.string()
                                             : event_args.columns_tsv_path;

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
    if (!output_path.parent_path().empty())
    {
        header.event_output_dir = output_path.parent_path().string();
    }
    if (!output_path.parent_path().empty())
    {
        std::filesystem::create_directories(output_path.parent_path());
    }

    const EventColumnProvider column_provider(
        {},
        {},
        columns_tsv_path);

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

    log_event_finish(log_prefix, entries.size(), elapsed_seconds);

    return 0;
}
