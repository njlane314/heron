/* -- C++ -- */
#ifndef NUXSEC_APPS_APP_UTILS_H
#define NUXSEC_APPS_APP_UTILS_H

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AppLog.hh"

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

inline const char *getenv_cstr(const char *name)
{
    const char *value = std::getenv(name);
    if (!value || !*value)
    {
        return nullptr;
    }
    return value;
}

inline std::filesystem::path repo_root_dir()
{
    if (const char *value = getenv_cstr("NUXSEC_REPO_ROOT"))
    {
        return std::filesystem::path(value);
    }
    return std::filesystem::current_path();
}

inline std::filesystem::path out_base_dir()
{
    if (const char *value = getenv_cstr("NUXSEC_OUT_BASE"))
    {
        return std::filesystem::path(value);
    }
    return repo_root_dir() / "scratch" / "out";
}

inline std::string workspace_set()
{
    if (const char *value = getenv_cstr("NUXSEC_SET"))
    {
        return std::string(value);
    }
    return "template";
}

inline std::filesystem::path stage_output_dir(const char *override_env, const std::string &stage)
{
    if (const char *value = getenv_cstr(override_env))
    {
        return std::filesystem::path(value);
    }
    return out_base_dir() / workspace_set() / stage;
}

inline int run_guarded(const std::string &log_prefix, const std::function<int()> &func)
{
    try
    {
        return func();
    }
    catch (const std::exception &e)
    {
        nuxsec::app::log::log_error(log_prefix, std::string("fatal_error=") + e.what());
        return 1;
    }
}

inline int run_guarded(const std::function<int()> &func)
{
    return run_guarded("nuxsec", func);
}

inline std::vector<std::string> read_paths(const std::string &filelist_path)
{
    std::ifstream fin(filelist_path);
    if (!fin)
    {
        throw std::runtime_error("Failed to open filelist: " + filelist_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) +
                                 "). Ensure the filelist exists (e.g. run scripts/partition-lists.sh).");
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

}

}

#endif
