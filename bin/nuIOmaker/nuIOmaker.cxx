/**
 *  @file  bin/nuIOmaker/nuIOmaker.cxx
 *
 *  @brief Main entrypoint for NuIO maker
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

#include "NuIO/NuIOMaker.h"

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

std::vector<std::string> read_file_list(const std::string &filelist_path)
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

struct Args
{
    std::string sample_name;
    std::string filelist_path;
    std::string output_path;
};

Args parse_args(int argc, char **argv)
{
    if (argc != 3)
    {
        throw std::runtime_error("Usage: nuIOmaker SAMPLE:FILELIST OUTPUT.root");
    }

    const std::string spec = argv[1];
    const auto pos = spec.find(':');
    if (pos == std::string::npos)
    {
        throw std::runtime_error("Bad sample spec (expected SAMPLE:FILELIST): " + spec);
    }

    Args args;
    args.sample_name = trim(spec.substr(0, pos));
    args.filelist_path = trim(spec.substr(pos + 1));
    args.output_path = argv[2];

    if (args.sample_name.empty() || args.filelist_path.empty() || args.output_path.empty())
    {
        throw std::runtime_error("Bad sample spec: " + spec);
    }

    return args;
}

}

int main(int argc, char **argv)
{
    using namespace nuio;

    try
    {
        const Args args = parse_args(argc, argv);
        const auto files = read_file_list(args.filelist_path);

        NuIOMakerConfig config;
        config.sample_name = args.sample_name;
        config.artio_files = files;
        config.output_path = args.output_path;

        NuIOMaker::write(config);

        std::cerr << "[nuIOmaker] wrote=" << config.output_path
                  << " sample=" << config.sample_name
                  << " files=" << config.artio_files.size()
                  << "\n";

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
