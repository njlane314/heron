/* -- C++ -- */
#ifndef NUXSEC_APPS_APP_UTILS_H
#define NUXSEC_APPS_APP_UTILS_H

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
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

inline int run_guarded(const std::function<int()> &func)
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

inline std::vector<std::string> read_paths(const std::string &filelist_path)
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

}

}

#endif
