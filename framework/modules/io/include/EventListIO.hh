/* -- C++ -- */
/**
 *  @file  framework/io/include/EventListIO.hh
 *
 *  @brief Event-level IO for selection and analysis bookkeeping outputs,
 *         capturing per-event metadata and processing summaries.
 */

#ifndef HERON_IO_EVENT_LIST_IO_H
#define HERON_IO_EVENT_LIST_IO_H

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "SampleIO.hh"

namespace nu
{
struct EventListHeader
{
    std::string analysis_name;
    std::string provenance_tree;
    std::string event_tree;
    std::string sample_list_source;
    std::string heron_set;
    std::string event_output_dir;
};

struct SampleInfo
{
    std::string sample_name;
    std::string sample_rootio_path;
    int sample_origin = -1;
    int beam_mode = -1;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;
};

class EventListIO
{
  public:
    enum class OpenMode
    {
        kRead,
        kUpdate
    };

    static void init(const std::string &out_path,
                     const EventListHeader &header,
                     const std::vector<SampleInfo> &sample_refs,
                     const std::string &event_schema_tsv,
                     const std::string &schema_tag);

    static EventListIO read(std::string path);

    explicit EventListIO(std::string path, OpenMode mode = OpenMode::kRead);

    const std::string &path() const noexcept { return m_path; }
    const EventListHeader &header() const noexcept { return m_header; }

    const std::unordered_map<int, SampleInfo> &sample_refs() const noexcept { return m_sample_refs; }

    std::string event_tree() const;

    ROOT::RDataFrame rdf() const;

    std::shared_ptr<const std::vector<char>> mask_for_origin(SampleIO::SampleOrigin origin) const;
    std::shared_ptr<const std::vector<char>> mask_for_mc_like() const;
    std::shared_ptr<const std::vector<char>> mask_for_data() const;
    std::shared_ptr<const std::vector<char>> mask_for_ext() const;

    double total_pot_data() const;
    double total_pot_mc() const;

    std::string beamline_label() const;

    std::string sample_tree_name(const std::string &sample_name,
                                 const std::string &tree_prefix = "events") const;

    ULong64_t snapshot_event_list(ROOT::RDF::RNode node,
                                  const std::string &sample_name,
                                  const std::vector<std::string> &columns,
                                  const std::string &selection = "true",
                                  const std::string &tree_prefix = "events",
                                  bool overwrite_if_exists = true) const;

    ULong64_t snapshot_event_list_merged(ROOT::RDF::RNode node,
                                         int sample_id,
                                         const std::string &sample_name,
                                         const std::vector<std::string> &columns,
                                         const std::string &selection,
                                         const std::string &tree_name = "events") const;

  private:
    std::string m_path;
    OpenMode m_mode;
    EventListHeader m_header{};
    std::unordered_map<int, SampleInfo> m_sample_refs;
    int m_max_sample_id = -1;
};
}

#endif
