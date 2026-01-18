/* -- C++ -- */
/**
 *  @file  apps/src/sampleIOaggregator.cc
 *
 *  @brief Main entrypoint for SampleIO provenance aggregation.
 */

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <tuple>
#include <stdexcept>
#include <string>
#include <vector>

#include "AppUtils.hh"
#include "SampleIO.hh"

namespace
{

struct Args
{
    std::string sample_name;
    std::string filelist_path;
    std::string output_path;
    std::string sample_list_path;
};

Args parse_args(int argc, char **argv)
{
    if (argc != 2)
    {
        throw std::runtime_error("Usage: sampleIOaggregator NAME:FILELIST");
    }

    const std::string spec = argv[1];
    const auto pos = spec.find(':');
    if (pos == std::string::npos)
    {
        throw std::runtime_error("Bad sample spec (expected NAME:FILELIST): " + spec);
    }

    Args args;
    args.sample_name = nuxsec::app::trim(spec.substr(0, pos));
    args.filelist_path = nuxsec::app::trim(spec.substr(pos + 1));

    if (args.sample_name.empty() || args.filelist_path.empty())
    {
        throw std::runtime_error("Bad sample spec: " + spec);
    }

    args.output_path = "./SampleIO_" + args.sample_name + ".root";
    args.sample_list_path = "./SampleIO_samples.tsv";

    return args;
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
        if (errno == ENOENT)
        {
            return {};
        }
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

void write_sample_list(const std::string &list_path, std::vector<SampleListEntry> entries)
{
    std::sort(entries.begin(), entries.end(),
              [](const SampleListEntry &a, const SampleListEntry &b)
              {
                  return std::tie(a.sample_kind, a.beam_mode, a.sample_name) <
                         std::tie(b.sample_kind, b.beam_mode, b.sample_name);
              });

    std::ofstream fout(list_path, std::ios::trunc);
    if (!fout)
    {
        throw std::runtime_error("Failed to open sample list for writing: " + list_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    fout << "# sample_name\tsample_kind\tbeam_mode\toutput_path\n";
    for (const auto &entry : entries)
    {
        fout << entry.sample_name << "\t"
             << entry.sample_kind << "\t"
             << entry.beam_mode << "\t"
             << entry.output_path << "\n";
    }
}

void update_sample_list(const std::string &list_path, const nuxsec::Sample &sample,
                        const std::string &output_path)
{
    auto entries = read_sample_list(list_path);
    const std::string kind_name = nuxsec::sample_kind_name(sample.kind);
    const std::string beam_name = nuxsec::beam_mode_name(sample.beam);

    bool updated = false;
    for (auto &entry : entries)
    {
        if (entry.sample_name == sample.sample_name &&
            entry.sample_kind == kind_name &&
            entry.beam_mode == beam_name)
        {
            entry.output_path = output_path;
            updated = true;
            break;
        }
    }
    if (!updated)
    {
        SampleListEntry entry;
        entry.sample_name = sample.sample_name;
        entry.sample_kind = kind_name;
        entry.beam_mode = beam_name;
        entry.output_path = output_path;
        entries.push_back(std::move(entry));
    }

    write_sample_list(list_path, std::move(entries));
}

}

int main(int argc, char **argv)
{
    using namespace nuxsec;

    try
    {
        const Args args = parse_args(argc, argv);
        const auto files = nuxsec::app::read_file_list(args.filelist_path);

        Sample sample = SampleIO::aggregate(args.sample_name, files);
        SampleIO::write(sample, args.output_path);
        update_sample_list(args.sample_list_path, sample, args.output_path);

        std::cerr << "[sampleIOaggregator] sample=" << sample.sample_name
                  << " stages=" << sample.stages.size()
                  << " pot_sum=" << sample.subrun_pot_sum
                  << " db_tortgt_pot_sum=" << sample.db_tortgt_pot_sum
                  << " normalization=" << sample.normalization
                  << " normalized_pot_sum=" << sample.normalized_pot_sum
                  << " output=" << args.output_path
                  << " sample_list=" << args.sample_list_path
                  << "\n";

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
