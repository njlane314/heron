/* -- C++ -- */
/**
 *  @file  framework/core/include/EventCLI.hh
 *
 *  @brief CLI helpers that drive event-level workflows, including configuration,
 *         input selection, and summary output for analysis-ready datasets.
 */
#ifndef HERON_CORE_EVENTCLI_H
#define HERON_CORE_EVENTCLI_H

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "AppLog.hh"
#include "SampleCLI.hh"
#include "SampleIO.hh"




inline void log_event_start(const std::string &log_prefix, const size_t sample_count)
{
    log_info(
        log_prefix,
        "action=event_build status=start samples=" +
            format_count(static_cast<long long>(sample_count)));
}

inline void log_event_finish(const std::string &log_prefix,
                             const size_t sample_count,
                             const double elapsed_seconds)
{
    std::ostringstream out;
    out << "action=event_build status=complete samples="
        << format_count(static_cast<long long>(sample_count))
        << " elapsed_s=" << std::fixed << std::setprecision(1)
        << elapsed_seconds;
    log_success(log_prefix, out.str());
}

void ensure_tree_present(const SampleIO::Sample &sample,
                         const std::string &tree_name);

struct EventArgs
{
    std::string list_path;
    std::string output_root;
    std::string selection;
    std::string columns_tsv_path;
};

struct EventInput
{
    SampleListEntry entry;
    SampleIO::Sample sample;
};

EventArgs parse_event_args(const std::vector<std::string> &args, const std::string &usage);

int run(const EventArgs &event_args, const std::string &log_prefix);




#endif // HERON_CORE_EVENTCLI_H
