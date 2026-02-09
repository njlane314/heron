/* -- C++ -- */
/**
 *  @file  io/include/EventList.hh
 *
 *  @brief Event list helper for filtering and metadata access.
 */

#ifndef NUXSEC_IO_EVENT_LIST_H
#define NUXSEC_IO_EVENT_LIST_H

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "EventIO.hh"   // Header, SampleInfo
#include "SampleIO.hh"  // SampleIO::SampleOrigin, BeamMode


class EventList
{
  public:
    explicit EventList(std::string path);

    const std::string &path() const noexcept { return m_path; }
    const Header &header() const noexcept { return m_header; }

    /// Sample refs, indexed by sample_id.
    const std::unordered_map<int, SampleInfo> &sample_refs() const noexcept { return m_sample_refs; }

    /// Tree name convenience (falls back to "events").
    std::string event_tree() const;

    /// Make an RDF on the event list tree.
    ROOT::RDataFrame rdf() const;

    /// Masks for filtering by sample_id (char is compact + trivially copyable).
    std::shared_ptr<const std::vector<char>> mask_for_origin(SampleIO::SampleOrigin origin) const;
    std::shared_ptr<const std::vector<char>> mask_for_mc_like() const;   // overlay + dirt + strangeness (+unknown if desired)
    std::shared_ptr<const std::vector<char>> mask_for_data() const;      // data only
    std::shared_ptr<const std::vector<char>> mask_for_ext() const;       // ext only

    /// Convenience POT totals (for watermark).
    double total_pot_data() const;
    double total_pot_mc() const;

    /// Beamline guess ("numi", "bnb", or "mixed"/"unknown").
    std::string beamline_label() const;

  private:
    std::string m_path;
    Header m_header{};
    std::unordered_map<int, SampleInfo> m_sample_refs;
    int m_max_sample_id = -1;
};

#endif
