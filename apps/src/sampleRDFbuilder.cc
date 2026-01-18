/* -- C++ -- */
/**
 *  @file  apps/src/sampleRDFbuilder.cc
 *
 *  @brief Build ROOT RDataFrame snapshots from SampleIO sample lists.
 */

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AppUtils.hh"
#include "AnalysisProcessor.hh"
#include "RDFBuilder.hh"
#include "SampleIO.hh"

namespace
{

struct SampleListEntry
{
    std::string sample_name;
    std::string sample_kind;
    std::string beam_mode;
    std::string output_path;
};

std::vector<SampleListEntry> read_sample_list(const std::string &list_path)
{
    std::ifstream fin(list_path);
    if (!fin)
    {
        throw std::runtime_error("Failed to open sample list: " + list_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }

    std::vector<SampleListEntry> entries;
    std::string line;
    while (std::getline(fin, line))
    {
        line = nuxsec::app::trim(line);
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        const auto fields = nuxsec::app::split_tabs(line);
        if (fields.size() < 4)
        {
            throw std::runtime_error("Malformed sample list entry: " + line);
        }
        SampleListEntry entry;
        entry.sample_name = fields[0];
        entry.sample_kind = fields[1];
        entry.beam_mode = fields[2];
        entry.output_path = fields[3];
        entries.push_back(std::move(entry));
    }
    return entries;
}

struct Args
{
    std::string list_path;
    std::string tree_name;
    std::string output_dir;
};

Args parse_args(int argc, char **argv)
{
    if (argc < 3 || argc > 4)
    {
        throw std::runtime_error("Usage: sampleRDFbuilder SAMPLE_LIST.tsv TREE_NAME [OUTPUT_DIR]");
    }

    Args args;
    args.list_path = nuxsec::app::trim(argv[1]);
    args.tree_name = nuxsec::app::trim(argv[2]);
    args.output_dir = (argc == 4) ? nuxsec::app::trim(argv[3]) : ".";

    if (args.list_path.empty() || args.tree_name.empty() || args.output_dir.empty())
    {
        throw std::runtime_error("Invalid arguments (empty path or tree name)");
    }

    return args;
}

std::string make_output_path(const std::string &output_dir, const std::string &sample_name)
{
    const std::string filename = "RDF_" + sample_name + ".root";
    std::filesystem::path out_dir(output_dir);
    return (out_dir / filename).string();
}

}

int main(int argc, char **argv)
{
    try
    {
        const Args args = parse_args(argc, argv);
        std::filesystem::create_directories(args.output_dir);

        const auto entries = read_sample_list(args.list_path);
        if (entries.empty())
        {
            throw std::runtime_error("Sample list is empty: " + args.list_path);
        }

        for (const auto &entry : entries)
        {
            const nuxsec::Sample sample = nuxsec::SampleIO::read(entry.output_path);
            ROOT::RDataFrame rdf = nuxsec::RDFBuilder::load_sample(sample, args.tree_name);
            nuxsec::ProcessorEntry proc_entry;
            switch (sample.kind)
            {
            case nuxsec::SampleKind::kData:
                proc_entry.source = nuxsec::SourceKind::kData;
                break;
            case nuxsec::SampleKind::kEXT:
                proc_entry.source = nuxsec::SourceKind::kExt;
                proc_entry.trig_nom = sample.db_tor101_pot_sum;
                proc_entry.trig_eqv = sample.subrun_pot_sum;
                break;
            case nuxsec::SampleKind::kMCOverlay:
            case nuxsec::SampleKind::kMCDirt:
            case nuxsec::SampleKind::kMCStrangeness:
                proc_entry.source = nuxsec::SourceKind::kMC;
                proc_entry.pot_nom = sample.db_tortgt_pot_sum;
                proc_entry.pot_eqv = sample.subrun_pot_sum;
                break;
            default:
                proc_entry.source = nuxsec::SourceKind::kUnknown;
                break;
            }

            const auto &processor = nuxsec::AnalysisProcessor::Processor();
            ROOT::RDF::RNode updated = processor.Run(rdf, proc_entry);
            const std::string output_path = make_output_path(args.output_dir, sample.sample_name);
            auto written = updated.Snapshot(args.tree_name, output_path);

            std::cerr << "[sampleRDFbuilder] sample=" << sample.sample_name
                      << " kind=" << entry.sample_kind
                      << " beam=" << entry.beam_mode
                      << " entries=" << written.GetValue()
                      << " output=" << output_path
                      << "\n";
        }

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
