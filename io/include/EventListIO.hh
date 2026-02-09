#ifndef NUXSEC_IO_EVENT_LIST_IO_H
#define NUXSEC_IO_EVENT_LIST_IO_H

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "EventIO.hh"
#include "SampleIO.hh"

namespace nu
{
class EventListIO
{
  public:
    explicit EventListIO(std::string path);

    const std::string &path() const noexcept { return m_path; }
    const Header &header() const noexcept { return m_header; }

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

    static int build_event_list(const std::string &out_root,
                                const std::string &samples_tsv,
                                const std::string &base_sel,
                                const std::vector<std::string> &columns);

  private:
    std::string m_path;
    Header m_header{};
    std::unordered_map<int, SampleInfo> m_sample_refs;
    int m_max_sample_id = -1;
};
}

#endif
