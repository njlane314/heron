/* -- C++ -- */
/**
 *  @file  io/include/EventIO.hh
 *
 *  @brief Event-level IO for selection/analysis bookkeeping outputs.
 */

#ifndef NUXSEC_IO_EVENTIO_H
#define NUXSEC_IO_EVENTIO_H

#include <cstdint>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

namespace nuxsec
{
namespace io
{

struct EventHeader
{
    std::string analysis_name;
    std::string analysis_tree;
    std::string sample_list_source;
};

struct EventSampleRef
{
    std::string sample_name;
    std::string sample_rootio_path;
    int sample_kind = -1;
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
                     const EventHeader &header,
                     const std::vector<EventSampleRef> &sample_refs,
                     const std::string &event_schema_tsv,
                     const std::string &schema_tag);

    explicit EventIO(std::string out_path, OpenMode mode);
    const std::string &path() const noexcept { return m_path; }

    std::string tree_name_for_sample(const std::string &sample_name,
                                     const std::string &tree_prefix = "events") const;

    ULong64_t snapshot_event_list(ROOT::RDF::RNode node,
                                  const std::string &sample_name,
                                  const std::string &x_column,
                                  const std::string &y_column,
                                  const std::vector<std::string> &double_columns,
                                  const std::vector<std::string> &int_columns,
                                  const std::string &selection = "true",
                                  const std::string &tree_prefix = "events",
                                  bool overwrite_if_exists = true) const;

private:
    static std::string sanitise_root_key(std::string s);

private:
    std::string m_path;
    OpenMode m_mode;
};

} // namespace io
} // namespace nuxsec

#endif // NUXSEC_IO_EVENTIO_H
