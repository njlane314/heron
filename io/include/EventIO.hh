/* -- C++ -- */
/**
 *  @file  io/include/EventIO.hh
 *
 *  @brief Event-level IO for selection/analysis bookkeeping outputs.
 */

#ifndef NUXSEC_EVENT_EVENTIO_H
#define NUXSEC_EVENT_EVENTIO_H

#include <cstdint>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>


struct Header
{
    std::string analysis_name;
    std::string provenance_tree;
    std::string event_tree;
    std::string sample_list_source;
    std::string nuxsec_set;
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

class EventIO final
{
public:
    enum class OpenMode
    {
        kRead,
        kUpdate
    };

    static void init(const std::string &out_path,
                     const Header &header,
                     const std::vector<SampleInfo> &sample_refs,
                     const std::string &event_schema_tsv,
                     const std::string &schema_tag);

    explicit EventIO(std::string out_path, OpenMode mode);
    const std::string &path() const noexcept { return m_path; }

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
};

}
}

#endif
