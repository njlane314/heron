/* -- C++ -- */
/**
 *  @file  framework/core/src/SampleCLI.cc
 *
 *  @brief Non-inline SampleCLI helpers for parsing and sample list I/O.
 */

#include "SampleCLI.hh"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "AppUtils.hh"

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

std::vector<SampleListEntry> read_samples(const std::string &list_path,
                                          bool allow_missing,
                                          bool require_nonempty)
{
    std::ifstream fin(list_path);
    if (!fin)
    {
        if (allow_missing && errno == ENOENT)
        {
            return {};
        }
        throw std::runtime_error("Failed to open sample list: " + list_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }

    std::vector<SampleListEntry> entries;
    std::string line;
    bool first_nonempty = true;
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

        if (first_nonempty && fields[0] == "sample_name")
        {
            first_nonempty = false;
            continue;
        }
        first_nonempty = false;

        SampleListEntry entry;
        entry.sample_name = fields[0];
        entry.sample_origin = fields[1];
        entry.beam_mode = fields[2];
        entry.output_path = fields[3];
        entries.push_back(std::move(entry));
    }

    if (require_nonempty && entries.empty())
    {
        throw std::runtime_error("Sample list is empty: " + list_path);
    }

    return entries;
}

void write_samples(const std::string &list_path, std::vector<SampleListEntry> entries)
{
    std::sort(entries.begin(), entries.end(),
              [](const SampleListEntry &a, const SampleListEntry &b)
              {
                  return std::tie(a.sample_origin, a.beam_mode, a.sample_name) <
                         std::tie(b.sample_origin, b.beam_mode, b.sample_name);
              });

    std::ofstream fout(list_path, std::ios::trunc);
    if (!fout)
    {
        throw std::runtime_error("Failed to open sample list for writing: " + list_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    fout << "# sample_name\tsample_origin\tbeam_mode\toutput_path\n";
    for (const auto &entry : entries)
    {
        fout << entry.sample_name << "\t"
             << entry.sample_origin << "\t"
             << entry.beam_mode << "\t"
             << entry.output_path << "\n";
    }
}

SampleArgs parse_sample_input(const std::string &input)
{
    const auto pos = input.find(':');
    if (pos == std::string::npos)
    {
        throw std::runtime_error("Bad sample definition (expected NAME:FILELIST): " + input);
    }

    SampleArgs out;
    out.sample_name = trim(input.substr(0, pos));
    out.filelist_path = trim(input.substr(pos + 1));

    if (out.sample_name.empty() || out.filelist_path.empty())
    {
        throw std::runtime_error("Bad sample definition: " + input);
    }

    const std::filesystem::path sample_dir =
        stage_output_dir("HERON_SAMPLE_DIR", "sample");
    out.output_path = (sample_dir / ("sample_root_" + out.sample_name + ".root")).string();
    out.sample_list_path = (sample_dir / "samples.tsv").string();

    return out;
}

SampleArgs parse_sample_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() != 1)
    {
        throw std::runtime_error(usage);
    }

    return parse_sample_input(args[0]);
}

void update_sample_list(const std::string &list_path,
                        const SampleIO::Sample &sample,
                        const std::string &output_path)
{
    auto entries = read_samples(list_path, true, false);
    const std::string origin_name = SampleIO::sample_origin_name(sample.origin);
    const std::string beam_name = SampleIO::beam_mode_name(sample.beam);

    bool updated = false;
    for (auto &entry : entries)
    {
        if (entry.sample_name == sample.sample_name &&
            entry.sample_origin == origin_name &&
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
        entry.sample_origin = origin_name;
        entry.beam_mode = beam_name;
        entry.output_path = output_path;
        entries.push_back(std::move(entry));
    }

    write_samples(list_path, std::move(entries));
}
