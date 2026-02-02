/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecEventIOdriver.cc
 *
 *  @brief Build event-level output from aggregated samples.
 */

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "AppUtils.hh"
#include "EventCLI.hh"
#include "EventColumnProvider.hh"
#include "RDataFrameService.hh"

namespace nuxsec
{

namespace app
{

int run(const event::Args &event_args, const std::string &log_prefix)
{
    ROOT::EnableImplicitMT();

    const auto &analysis = nuxsec::AnalysisConfigService::instance();
    const auto entries = nuxsec::app::sample::read_samples(event_args.list_path);
    
    const auto start_time = std::chrono::steady_clock::now();
    nuxsec::app::event::log_event_start(log_prefix, entries.size());

    nuxsec::app::StatusMonitor status_monitor(
        log_prefix,
        "action=event_build status=running message=processing");

    std::vector<nuxsec::app::event::Input> inputs;
    inputs.reserve(entries.size());
    
    std::vector<nuxsec::event::SampleInfo> sample_infos;
    sample_infos.reserve(entries.size());

    for (const auto &entry : entries)
    {
        nuxsec::sample::SampleIO::Sample sample = nuxsec::sample::SampleIO::read(entry.output_path);

        nuxsec::event::SampleInfo info;
        info.sample_name = sample.sample_name;
        info.sample_rootio_path = entry.output_path;
        info.sample_origin = static_cast<int>(sample.origin);
        info.beam_mode = static_cast<int>(sample.beam);
        info.subrun_pot_sum = sample.subrun_pot_sum;
        info.db_tortgt_pot_sum = sample.db_tortgt_pot_sum;
        info.db_tor101_pot_sum = sample.db_tor101_pot_sum;
        sample_infos.push_back(std::move(info));

        nuxsec::app::event::Input input;
        input.entry = entry;
        input.sample = std::move(sample);
        inputs.push_back(std::move(input));
    }

    const std::vector<std::string> default_columns = {
        "run",
        "sub",
        "evt",
        "is_signal",
        "w_nominal",
        "detector_image_u",
        "detector_image_v",
        "detector_image_w"
    };

    const std::vector<std::pair<std::string, std::string>> default_schema_columns = {
        {"int", "run"},
        {"int", "sub"},
        {"int", "evt"},
        {"bool", "is_signal"},
        {"double", "w_nominal"},
        {"vector<float>", "detector_image_u"},
        {"vector<float>", "detector_image_v"},
        {"vector<float>", "detector_image_w"}
    };

    const std::string provenance_tree = "nuxsec_art_provenance/run_subrun";
    const std::string event_tree = analysis.tree_name();

    nuxsec::event::Header header;
    header.analysis_name = analysis.name();
    header.provenance_tree = provenance_tree;
    header.event_tree = event_tree;
    header.sample_list_source = event_args.list_path;
    header.nuxsec_set = nuxsec::app::workspace_set();

    const std::filesystem::path output_path(event_args.output_root);
    if (!output_path.parent_path().empty())
    {
        header.event_output_dir = output_path.parent_path().string();
    }
    if (!output_path.parent_path().empty())
        std::filesystem::create_directories(output_path.parent_path());

    const nuxsec::app::event::EventColumnProvider column_provider(
        default_columns,
        default_schema_columns,
        event_args.columns_tsv_path);

    nuxsec::event::EventIO::init(event_args.output_root,
                                 header,
                                 sample_infos,
                                 column_provider.schema_tsv(),
                                 column_provider.schema_tag());
    nuxsec::event::EventIO event_io(event_args.output_root,
                                    nuxsec::event::EventIO::OpenMode::kUpdate);

    for (size_t i = 0; i < inputs.size(); ++i)
    {
        const auto &input = inputs[i];
        const nuxsec::sample::SampleIO::Sample &sample = input.sample;
        const int sample_id = static_cast<int>(i);
        
        nuxsec::app::log::log_stage(
            log_prefix,
            "ensure_tree",
            "sample=" + sample.sample_name + " tree=" + event_tree);

        nuxsec::app::event::ensure_tree_present(sample, event_tree);

        nuxsec::app::log::log_stage(
            log_prefix,
            "load_rdf",
            "sample=" + sample.sample_name);

        ROOT::RDataFrame rdf = nuxsec::RDataFrameService::load_sample(sample, event_tree);

        nuxsec::app::log::log_stage(
            log_prefix,
            "make_processor",
            "sample=" + sample.sample_name);

        const nuxsec::ProcessorEntry proc_entry = analysis.make_processor(sample);
        const auto &processor = nuxsec::ColumnDerivationService::instance();

        nuxsec::app::log::log_stage(
            log_prefix,
            "define_columns",
            "sample=" + sample.sample_name);

        ROOT::RDF::RNode node = processor.define(rdf, proc_entry);

        using SampleOrigin = nuxsec::sample::SampleIO::SampleOrigin;
        const auto origin = sample.origin;
        
        if (origin == SampleOrigin::kOverlay)
        {
            nuxsec::app::log::log_stage(
                log_prefix,
                "filter_overlay",
                "sample=" + sample.sample_name);
            node = node.Filter([](int strange) { return strange == 0; }, {"count_strange"});
        }
        else if (origin == SampleOrigin::kStrangeness)
        {
            nuxsec::app::log::log_stage(
                log_prefix,
                "filter_strangeness",
                "sample=" + sample.sample_name);
            node = node.Filter([](int strange) { return strange > 0; }, {"count_strange"});
        }

        std::string snapshot_message = "sample=" + sample.sample_name;
        if (!event_args.selection.empty())
        {
            snapshot_message += " selection=" + event_args.selection;
        }
        nuxsec::app::log::log_stage(
            log_prefix,
            "snapshot",
            snapshot_message);

        const ULong64_t n_written =
            event_io.snapshot_event_list_merged(node,
                                                sample_id,
                                                sample.sample_name,
                                                column_provider.columns(),
                                                event_args.selection,
                                                "events");

        std::ostringstream log_message;
        log_message << "action=event_snapshot status=complete analysis=" << analysis.name()
                    << " sample=" << sample.sample_name
                    << " kind=" << nuxsec::sample::SampleIO::sample_origin_name(sample.origin)
                    << " beam=" << nuxsec::sample::SampleIO::beam_mode_name(sample.beam)
                    << " events_written=" << n_written
                    << " output=" << event_args.output_root;
        if (!event_args.selection.empty())
        {
            log_message << " selection=" << event_args.selection;
        }
        nuxsec::app::log::log_success(log_prefix, log_message.str());
    }
    status_monitor.stop();

    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    
    nuxsec::app::event::log_event_finish(log_prefix, entries.size(), elapsed_seconds);

    return 0;
}

} // namespace app

} // namespace nuxsec

int main(int argc, char **argv)
{
    return nuxsec::app::run_guarded(
        "nuxsecEventIOdriver",
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::event::Args event_args =
                nuxsec::app::event::parse_args(
                    args,
                    "Usage: nuxsecEventIOdriver SAMPLE_LIST.tsv OUTPUT.root [SELECTION] [COLUMNS.tsv]");
            return nuxsec::app::run(event_args, "nuxsecEventIOdriver");
        });
}
