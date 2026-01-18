/**
 *  @file  bin/sampleIOaggregator/sampleIOaggregator.cxx
 *
 *  @brief Main entrypoint for SampleIO provenance aggregation
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
    args.sample_name = trim(spec.substr(0, pos));
    args.filelist_path = trim(spec.substr(pos + 1));

    if (args.sample_name.empty() || args.filelist_path.empty())
    {
        throw std::runtime_error("Bad sample spec: " + spec);
    }

    args.output_path = "./SampleIO_" + args.sample_name + ".root";

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

        Sample sample = SampleIO::aggregate(args.sample_name, files);
        SampleIO::write(sample, args.output_path);

        std::cerr << "[sampleIOaggregator] sample=" << sample.sample_name
                  << " stages=" << sample.stages.size()
                  << " pot_sum=" << sample.subrun_pot_sum
                  << " db_tortgt_pot_sum=" << sample.db_tortgt_pot_sum
                  << " normalization=" << sample.normalization
                  << " normalized_pot_sum=" << sample.normalized_pot_sum
                  << " output=" << args.output_path
                  << "\n";

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
