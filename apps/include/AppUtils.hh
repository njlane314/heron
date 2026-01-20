/* -- C++ -- */
#ifndef Nuxsec_APPS_APPUTILS_H_INCLUDED
#define Nuxsec_APPS_APPUTILS_H_INCLUDED

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace nuxsec
{

namespace app
{

inline std::string trim(std::string s)
{
    auto notspace = [](unsigned char c)
    {
        return std::isspace(c) == 0;
    };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
    return s;
}

inline std::vector<std::string> collect_args(int argc, char **argv, int start_index = 1)
{
    std::vector<std::string> args;
    if (argc <= start_index)
    {
        return args;
    }
    args.reserve(static_cast<size_t>(argc - start_index));
    for (int i = start_index; i < argc; ++i)
    {
        args.emplace_back(argv[i]);
    }
    return args;
}

inline int run_with_exceptions(const std::function<int()> &func)
{
    try
    {
        return func();
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}

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

inline std::vector<std::string> read_file_list(const std::string &filelist_path)
{
    std::ifstream fin(filelist_path);
    if (!fin)
    {
        throw std::runtime_error("Failed to open filelist: " + filelist_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    std::vector<std::string> files;
    std::string line;
    while (std::getline(fin, line))
    {
        line = trim(line);
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        files.push_back(line);
    }
    if (files.empty())
    {
        throw std::runtime_error("Filelist is empty: " + filelist_path);
    }
    return files;
}

struct SampleListEntry
{
    std::string sample_name;
    std::string sample_kind;
    std::string beam_mode;
    std::string output_path;
};

inline std::vector<SampleListEntry> read_sample_list(const std::string &list_path,
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

inline void write_sample_list(const std::string &list_path, std::vector<SampleListEntry> entries)
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

}

}

#endif
