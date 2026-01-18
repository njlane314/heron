/**
 *  @file  bin/sampleRDFmaker/sampleRDFmaker.cxx
 *
 *  @brief Main entrypoint for building ROOT RDataFrame from sample lists
 */

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <exception>
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

SampleListEntry find_sample_entry(const std::vector<SampleListEntry> &entries,
                                  const std::string &sample_name)
{
    std::vector<SampleListEntry> matches;
    for (const auto &entry : entries)
    {
        if (entry.sample_name == sample_name)
        {
            matches.push_back(entry);
        }
    }
    if (matches.empty())
    {
        throw std::runtime_error("Sample not found in list: " + sample_name);
    }
    if (matches.size() > 1)
    {
        throw std::runtime_error("Sample name is not unique in list: " + sample_name);
    }
    return matches.front();
}

struct Args
{
    std::string sample_list_path;
    std::string sample_name;
    std::string tree_name;
    std::string output_path;
    bool write_snapshot = false;
};

Args parse_args(int argc, char **argv)
{
    if (argc != 4 && argc != 5)
    {
        throw std::runtime_error(
            "Usage: sampleRDFmaker SAMPLE_LIST.tsv SAMPLE_NAME TREE_NAME [OUTPUT.root]");
    }

    Args args;
    args.sample_list_path = argv[1];
    args.sample_name = argv[2];
    args.tree_name = argv[3];

    if (argc == 5)
    {
        args.output_path = argv[4];
        args.write_snapshot = true;
    }

    if (args.sample_list_path.empty() || args.sample_name.empty() || args.tree_name.empty())
    {
        throw std::runtime_error("Sample list, sample name, and tree name are required");
    }

    return args;
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        const Args args = parse_args(argc, argv);
        const auto entries = read_sample_list(args.sample_list_path);
        const auto entry = find_sample_entry(entries, args.sample_name);

        const nuio::Sample sample = nuio::SampleIO::read(entry.output_path);
        ROOT::RDataFrame rdf = nuana::SampleRDF::load_sample(sample, args.tree_name);

        if (args.write_snapshot)
        {
            rdf.Snapshot(args.tree_name, args.output_path);
            std::cerr << "[sampleRDFmaker] wrote=" << args.output_path
                      << " sample=" << sample.sample_name
                      << " stages=" << sample.stages.size()
                      << " tree=" << args.tree_name
                      << "\n";
        }
        else
        {
            const auto count = rdf.Count();
            const auto total = count.GetValue();
            std::cerr << "[sampleRDFmaker] sample=" << sample.sample_name
                      << " stages=" << sample.stages.size()
                      << " tree=" << args.tree_name
                      << " entries=" << total
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
