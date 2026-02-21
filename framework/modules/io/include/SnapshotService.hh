/* -- C++ -- */
/**
 *  @file  framework/io/include/SnapshotService.hh
 *
 *  @brief Snapshot helpers for event-level ROOT output.
 */

#ifndef HERON_IO_SNAPSHOT_SERVICE_H
#define HERON_IO_SNAPSHOT_SERVICE_H

#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>


class SnapshotService final
{
  public:
    static std::string sanitise_root_key(std::string s);

    static ULong64_t snapshot_event_list(ROOT::RDF::RNode node,
                                         const std::string &out_path,
                                         const std::string &sample_name,
                                         const std::vector<std::string> &columns,
                                         const std::string &selection = "true",
                                         const std::string &tree_prefix = "events",
                                         bool overwrite_if_exists = true);

    static ULong64_t snapshot_event_list_merged(ROOT::RDF::RNode node,
                                                const std::string &out_path,
                                                int sample_id,
                                                const std::string &sample_name,
                                                const std::vector<std::string> &columns,
                                                const std::string &selection,
                                                const std::string &tree_name = "events");
};


#endif // HERON_IO_SNAPSHOT_SERVICE_H
