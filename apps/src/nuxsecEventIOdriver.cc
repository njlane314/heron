/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecEventIOdriver.cc
 *
 *  @brief Build event-level output from aggregated samples.
 */

#include "EventCLI.hh"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "RDataFrameService.hh"

namespace nuxsec
{

namespace app
{

int run(const event::Args &event_args, const std::string &log_prefix)
{
    ROOT::EnableImplicitMT();

    const auto &analysis = nuxsec::AnalysisConfigService::instance();
    const auto entries = nuxsec::app::read_samples(event_args.list_path);
    
    const auto start_time = std::chrono::steady_clock::now();
    nuxsec::app::event::log_event_start(log_prefix, entries.size());

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

    const std::vector<std::string> columns = {
        "run",
        "sub",
        "evt"
    };

    const std::string provenance_tree = "nuxsec_art_provenance/run_subrun";
    const std::string event_tree = analysis.tree_name();

    nuxsec::event::Header header;
    header.analysis_name = analysis.name();
    header.provenance_tree = provenance_tree;
    header.event_tree = event_tree;
    header.sample_list_source = event_args.list_path;

    const std::filesystem::path output_path(event_args.output_root);
    if (!output_path.parent_path().empty())
        std::filesystem::create_directories(output_path.parent_path());

    std::ostringstream schema;
    schema << "type\tname\n";
    for (const auto &c : columns)
        schema << "int\t" << c << "\n";

    nuxsec::event::EventIO::init(event_args.output_root,
                                 header,
                                 sample_infos,
                                 schema.str(),
                                 "compiled");
    nuxsec::event::EventIO event_io(event_args.output_root,
                                    nuxsec::event::EventIO::OpenMode::kUpdate);

    for (size_t i = 0; i < inputs.size(); ++i)
    {
        const auto &input = inputs[i];
        const nuxsec::sample::SampleIO::Sample &sample = input.sample;
        const int sample_id = static_cast<int>(i);
        
        std::cerr << "[" << log_prefix << "]"
                  << " stage=ensure_tree sample=" << sample.sample_name
                  << " tree=" << event_tree
                  << "\n";
        
        nuxsec::app::event::ensure_tree_present(sample, event_tree);
        
        std::cerr << "[" << log_prefix << "]"
                  << " stage=load_rdf sample=" << sample.sample_name
                  << "\n";
        
        ROOT::RDataFrame rdf = nuxsec::RDataFrameService::load_sample(sample, event_tree);
        
        std::cerr << "[" << log_prefix << "]"
                  << " stage=make_processor sample=" << sample.sample_name
                  << "\n";
        
        const nuxsec::ProcessorEntry proc_entry = analysis.make_processor(sample);
        const auto &processor = nuxsec::ColumnDerivationService::instance();
        
        std::cerr << "[" << log_prefix << "]"
                  << " stage=define_columns sample=" << sample.sample_name
                  << "\n";
        
        ROOT::RDF::RNode node = processor.define(rdf, proc_entry);

        using SampleOrigin = nuxsec::sample::SampleIO::SampleOrigin;
        const auto origin = sample.origin;
        
        if (origin == SampleOrigin::kOverlay)
        {
            std::cerr << "[" << log_prefix << "]"
                      << " stage=filter_overlay sample=" << sample.sample_name
                      << "\n";
            node = node.Filter([](int strange) { return strange == 0; }, {"count_strange"});
        }
        else if (origin == SampleOrigin::kStrangeness)
        {
            std::cerr << "[" << log_prefix << "]"
                      << " stage=filter_strangeness sample=" << sample.sample_name
                      << "\n";
            node = node.Filter([](int strange) { return strange > 0; }, {"count_strange"});
        }

        std::cerr << "[" << log_prefix << "]"
                  << " stage=snapshot sample=" << sample.sample_name
                  << "\n";

        const ULong64_t n_written =
            event_io.snapshot_event_list_merged(node,
                                                sample_id,
                                                sample.sample_name,
                                                columns,
                                                "",
                                                "events");

        std::cerr << "[" << log_prefix << "] analysis=" << analysis.name()
                  << " sample=" << sample.sample_name
                  << " kind=" << nuxsec::sample::SampleIO::sample_origin_name(sample.origin)
                  << " beam=" << nuxsec::sample::SampleIO::beam_mode_name(sample.beam)
                  << " events_written=" << n_written
                  << " output=" << event_args.output_root
                  << "\n";
    }
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
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::event::Args event_args =
                nuxsec::app::event::parse_args(
                    args, "Usage: nuxsecEventIOdriver SAMPLE_LIST.tsv OUTPUT.root");
            return nuxsec::app::run(event_args, "nuxsecEventIOdriver");
        });
}
