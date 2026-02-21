/* -- C++ -- */
/**
 *  @file  framework/core/include/SampleCLI.hh
 *
 *  @brief CLI helpers that manage sample-level workflows, from input handling
 *         through reporting and normalisation for data preparation tasks.
 */
#ifndef HERON_CORE_SAMPLECLI_H
#define HERON_CORE_SAMPLECLI_H

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "AppLog.hh"
#include "SampleIO.hh"

std::vector<std::string> split_tabs(const std::string &line);

struct SampleListEntry
{
    std::string sample_name;
    std::string sample_origin;
    std::string beam_mode;
    std::string output_path;
};

std::vector<SampleListEntry> read_samples(const std::string &list_path,
                                          bool allow_missing = false,
                                          bool require_nonempty = true);

void write_samples(const std::string &list_path, std::vector<SampleListEntry> entries);

inline void log_sample_start(const std::string &log_prefix, const size_t file_count)
{
    log_info(
        log_prefix,
        "action=sample_build status=start files=" +
            format_count(static_cast<long long>(file_count)));
}

inline void log_sample_finish(const std::string &log_prefix,
                              const size_t input_count,
                              const double elapsed_seconds)
{
    std::ostringstream out;
    out << "action=sample_build status=complete inputs="
        << format_count(static_cast<long long>(input_count))
        << " elapsed_s=" << std::fixed << std::setprecision(1)
        << elapsed_seconds;
    log_success(log_prefix, out.str());
}

struct SampleArgs
{
    std::string sample_name;
    std::string filelist_path;
    std::string output_path;
    std::string sample_list_path;
};

SampleArgs parse_sample_input(const std::string &input);

SampleArgs parse_sample_args(const std::vector<std::string> &args, const std::string &usage);

void update_sample_list(const std::string &list_path,
                        const SampleIO::Sample &sample,
                        const std::string &output_path);

int run(const SampleArgs &sample_args, const std::string &log_prefix);

#endif
