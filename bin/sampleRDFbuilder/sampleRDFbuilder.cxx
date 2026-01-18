/**
 *  @file  bin/sampleRDFbuilder/sampleRDFbuilder.cxx
 *
 *  @brief Build ROOT RDataFrame snapshots from SampleIO sample lists
 */

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "NuAna/SampleRDF.h"
#include "NuIO/SampleIO.h"

namespace
{

std::string trim(std::string s)
{
    auto notspace = [](unsigned char c)
    {
        return std::isspace(c) == 0;
    };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
    return s;
}

std::vector<std::string> split_tabs(const std::string &line)
{
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= line.size())
    {
        const size_t pos = line.find('\t', start);
        if (pos == std::string::npos)
        {
            out.push_back(line.substr(start));
            break;
        }
        out.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

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
        line = trim(line);
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        const auto fields = split_tabs(line);
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
    args.list_path = trim(argv[1]);
    args.tree_name = trim(argv[2]);
    args.output_dir = (argc == 4) ? trim(argv[3]) : ".";

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
            const nuio::Sample sample = nuio::SampleIO::read(entry.output_path);
            ROOT::RDataFrame rdf = nuana::SampleRDF::load_sample(sample, args.tree_name);
            const std::string output_path = make_output_path(args.output_dir, sample.sample_name);
            auto written = rdf.Snapshot(args.tree_name, output_path);

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
