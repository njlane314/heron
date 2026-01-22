/* -- C++ -- */
#ifndef NUXSEC_APPS_EVENTCOMMAND_H
#define NUXSEC_APPS_EVENTCOMMAND_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

#include "AnalysisConfigService.hh"
#include "ColumnDerivationService.hh"
#include "EventIO.hh"
#include "RDataFrameFactory.hh"
#include "SampleIO.hh"
#include "AppUtils.hh"

namespace nuxsec
{

namespace app
{

struct EventIOArgs
{
    std::string list_path;
    std::string output_root;
    int nthreads = 0;
};

struct EventIOInput
{
    nuxsec::app::SampleListEntry entry;
    sample::SampleIO::Sample sample;
};

inline EventIOArgs parse_eventio_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() != 2 && args.size() != 3)
    {
        throw std::runtime_error(usage);
    }

    EventIOArgs out;
    out.list_path = nuxsec::app::trim(args.at(0));
    out.output_root = nuxsec::app::trim(args.at(1));
    if (args.size() == 3)
    {
        try
        {
            out.nthreads = std::stoi(nuxsec::app::trim(args.at(2)));
        }
        catch (const std::exception &)
        {
            throw std::runtime_error("Invalid NTHREADS (expected integer)");
        }
        if (out.nthreads < 0)
            throw std::runtime_error("Invalid NTHREADS (must be >= 0)");
    }

    if (out.list_path.empty() || out.output_root.empty())
    {
        throw std::runtime_error("Invalid arguments (empty path)");
    }

    return out;
}

inline int run_eventio(const EventIOArgs &event_args, const std::string &log_prefix)
{
    if (event_args.nthreads > 0)
        ROOT::EnableImplicitMT(event_args.nthreads);
    else
        ROOT::EnableImplicitMT();

    const auto &analysis = nuxsec::AnalysisConfigService::instance();
    const auto entries = nuxsec::app::read_sample_list(event_args.list_path);

    std::vector<EventIOInput> inputs;
    inputs.reserve(entries.size());
    std::vector<nuxsec::io::EventSampleRef> sample_refs;
    sample_refs.reserve(entries.size());

    for (const auto &entry : entries)
    {
        sample::SampleIO::Sample sample = sample::SampleIO::read(entry.output_path);

        nuxsec::io::EventSampleRef ref;
        ref.sample_name = sample.sample_name;
        ref.sample_rootio_path = entry.output_path;
        ref.sample_kind = static_cast<int>(sample.kind);
        ref.beam_mode = static_cast<int>(sample.beam);
        ref.subrun_pot_sum = sample.subrun_pot_sum;
        ref.db_tortgt_pot_sum = sample.db_tortgt_pot_sum;
        ref.db_tor101_pot_sum = sample.db_tor101_pot_sum;
        sample_refs.push_back(std::move(ref));

        EventIOInput input;
        input.entry = entry;
        input.sample = std::move(sample);
        inputs.push_back(std::move(input));
    }

    const std::string x_column = "evt_x";
    const std::string y_column;
    const std::vector<std::string> double_columns = {
        "w_nominal",
        "w_template",
        "reco_nu_energy"
    };
    const std::vector<std::string> int_columns = {
        "run_i",
        "subrun_i",
        "event_i",
        "evt_cutmask_i",
        "evt_category",
        "sel_template_i",
        "sel_reco_fv_i",
        "sel_signal_i",
        "sel_bkg_i",
        "in_reco_fiducial_i",
        "analysis_channels",
        "scattering_mode",
        "is_strange_i",
        "count_strange",
        "is_signal_i",
        "recognised_signal_i"
    };

    std::ostringstream schema;
    schema << "type\tname\n";
    schema << "x\t" << x_column << "\n";
    if (!y_column.empty())
        schema << "y\t" << y_column << "\n";
    for (const auto &c : double_columns)
        schema << "double\t" << c << "\n";
    for (const auto &c : int_columns)
        schema << "int\t" << c << "\n";

    nuxsec::io::EventHeader header;
    header.analysis_name = analysis.name();
    header.analysis_tree = analysis.tree_name();
    header.sample_list_source = event_args.list_path;

    nuxsec::io::EventIO::init(event_args.output_root, header, sample_refs, schema.str(), "compiled");
    nuxsec::io::EventIO event_io(event_args.output_root, nuxsec::io::EventIO::OpenMode::kUpdate);

    for (const auto &input : inputs)
    {
        const sample::SampleIO::Sample &sample = input.sample;
        ROOT::RDataFrame rdf = nuxsec::RDataFrameFactory::load_sample(sample, analysis.tree_name());
        const nuxsec::ProcessorEntry proc_entry = analysis.make_processor_entry(sample);

        const auto &processor = nuxsec::ColumnDerivationService::instance();
        ROOT::RDF::RNode node = processor.define(rdf, proc_entry);

        {
            const auto cnames = node.GetColumnNames();
            auto has = [&](const std::string &name) {
                return std::find(cnames.begin(), cnames.end(), name) != cnames.end();
            };
            if (!has("run"))
                node = node.Define("run", [] { return static_cast<unsigned int>(0); });
            if (!has("subrun"))
                node = node.Define("subrun", [] { return static_cast<unsigned int>(0); });
            if (!has("event"))
                node = node.Define("event", [] { return static_cast<unsigned int>(0); });

            if (!has("reco_nu_energy"))
                node = node.Define("reco_nu_energy", [] { return 0.0f; });
        }

        node = node.Define(
            "evt_cutmask",
            [](bool t, bool fv, bool sig) {
                uint64_t m = 0;
                if (t)   m |= (1ull << 0);
                if (fv)  m |= (1ull << 1);
                if (sig) m |= (1ull << 2);
                return m;
            },
            {"sel_template", "in_reco_fiducial", "recognised_signal"});

        node = node.Define(
            "evt_category",
            [](bool s, bool b) { return s ? 1 : (b ? 0 : -1); },
            {"sel_signal", "sel_bkg"});

        node = node.Define("evt_x", [](float e) { return static_cast<double>(e); }, {"reco_nu_energy"});
        node = node.Define("run_i", [](unsigned int v) { return static_cast<int>(v); }, {"run"});
        node = node.Define("subrun_i", [](unsigned int v) { return static_cast<int>(v); }, {"subrun"});
        node = node.Define("event_i", [](unsigned int v) { return static_cast<int>(v); }, {"event"});
        node = node.Define("evt_cutmask_i", [](uint64_t v) { return static_cast<int>(v); }, {"evt_cutmask"});
        node = node.Define("sel_template_i", [](bool v) { return static_cast<int>(v); }, {"sel_template"});
        node = node.Define("sel_reco_fv_i", [](bool v) { return static_cast<int>(v); }, {"sel_reco_fv"});
        node = node.Define("sel_signal_i", [](bool v) { return static_cast<int>(v); }, {"sel_signal"});
        node = node.Define("sel_bkg_i", [](bool v) { return static_cast<int>(v); }, {"sel_bkg"});
        node = node.Define("in_reco_fiducial_i", [](bool v) { return static_cast<int>(v); }, {"in_reco_fiducial"});
        node = node.Define("is_strange_i", [](bool v) { return static_cast<int>(v); }, {"is_strange"});
        node = node.Define("is_signal_i", [](bool v) { return static_cast<int>(v); }, {"is_signal"});
        node = node.Define("recognised_signal_i", [](bool v) { return static_cast<int>(v); }, {"recognised_signal"});

        const ULong64_t n_written =
            event_io.snapshot_event_list(node,
                                         sample.sample_name,
                                         x_column,
                                         y_column,
                                         double_columns,
                                         int_columns,
                                         "sel_template");

        std::cerr << "[" << log_prefix << "] analysis=" << analysis.name()
                  << " sample=" << sample.sample_name
                  << " kind=" << sample::SampleIO::sample_kind_name(sample.kind)
                  << " beam=" << sample::SampleIO::beam_mode_name(sample.beam)
                  << " events_written=" << n_written
                  << " output=" << event_args.output_root
                  << "\n";
    }

    return 0;
}

} // namespace app

} // namespace nuxsec

#endif // NUXSEC_APPS_EVENTCOMMAND_H
