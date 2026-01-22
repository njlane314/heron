/* -- C++ -- */
#ifndef NUXSEC_APPS_SAMPLE_UTILS_H
#define NUXSEC_APPS_SAMPLE_UTILS_H

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AppUtils.hh"

namespace nuxsec
{

namespace app
{

inline std::vector<std::string> split_tabs(const std::string &line)
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

inline std::vector<SampleListEntry> read_samples(const std::string &list_path,
                                                 bool allow_missing = false,
                                                 bool require_nonempty = true)
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
        entry.sample_kind = fields[1];
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

inline void write_samples(const std::string &list_path, std::vector<SampleListEntry> entries)
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
    fout << "# sample_name\tsample_origin\tbeam_mode\toutput_path\n";
    for (const auto &entry : entries)
    {
        fout << entry.sample_name << "\t"
             << entry.sample_kind << "\t"
             << entry.beam_mode << "\t"
             << entry.output_path << "\n";
    }
}

}

}

#endif
