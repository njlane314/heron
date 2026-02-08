/* -- C++ -- */
/**
 *  @file  io/src/EventIO.cc
 *
 *  @brief Implementation of event-level IO.
 */

#include "EventIO.hh"
#include "SnapshotService.hh"

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <system_error>

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

void EventIO::init(const std::string &out_path,
                   const Header &header,
                   const std::vector<SampleInfo> &sample_refs,
                   const std::string &event_schema_tsv,
                   const std::string &schema_tag)
{
    const std::filesystem::path output_path(out_path);
    if (!output_path.parent_path().empty())
    {
        std::error_code ec;
        std::filesystem::create_directories(output_path.parent_path(), ec);
        if (ec)
        {
            throw std::runtime_error(
                "EventIO::init: failed to create output directory: "
                + output_path.parent_path().string()
                + " (" + ec.message() + ")");
        }
    }

    std::unique_ptr<TFile> fout(TFile::Open(out_path.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie())
        throw std::runtime_error("EventIO::init: failed to create output file: " + out_path);

    TObjString(header.analysis_name.c_str()).Write("analysis_name");
    TObjString(header.provenance_tree.c_str()).Write("provenance_tree");
    TObjString(header.event_tree.c_str()).Write("event_tree");
    TObjString(header.sample_list_source.c_str()).Write("sample_list_source");
    if (!header.nuxsec_set.empty())
    {
        TObjString(header.nuxsec_set.c_str()).Write("nuxsec_set");
    }
    if (!header.event_output_dir.empty())
    {
        TObjString(header.event_output_dir.c_str()).Write("event_output_dir");
    }

    if (!event_schema_tsv.empty())
    {
        const std::string key = schema_tag.empty()
                                    ? "event_schema"
                                    : ("event_schema_" + SnapshotService::sanitise_root_key(schema_tag));
        TObjString(event_schema_tsv.c_str()).Write(key.c_str());
    }

    TTree tref("sample_refs", "Sample references (source + POT/triggers)");
    int sample_id = -1;
    std::string sample_name;
    std::string sample_rootio_path;
    int sample_origin = -1;
    int beam_mode = -1;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;

    tref.Branch("sample_id", &sample_id);
    tref.Branch("sample_name", &sample_name);
    tref.Branch("sample_rootio_path", &sample_rootio_path);
    tref.Branch("sample_origin", &sample_origin);
    tref.Branch("beam_mode", &beam_mode);
    tref.Branch("subrun_pot_sum", &subrun_pot_sum);
    tref.Branch("db_tortgt_pot_sum", &db_tortgt_pot_sum);
    tref.Branch("db_tor101_pot_sum", &db_tor101_pot_sum);

    for (size_t i = 0; i < sample_refs.size(); ++i)
    {
        const auto &r = sample_refs[i];
        sample_id = static_cast<int>(i);
        sample_name = r.sample_name;
        sample_rootio_path = r.sample_rootio_path;
        sample_origin = r.sample_origin;
        beam_mode = r.beam_mode;
        subrun_pot_sum = r.subrun_pot_sum;
        db_tortgt_pot_sum = r.db_tortgt_pot_sum;
        db_tor101_pot_sum = r.db_tor101_pot_sum;
        tref.Fill();
    }

    tref.Write();
    fout->Close();
}

EventIO::EventIO(std::string out_path, OpenMode mode)
    : m_path(std::move(out_path)), m_mode(mode)
{
    const char *opt = (m_mode == OpenMode::kRead) ? "READ" : "UPDATE";
    std::unique_ptr<TFile> f(TFile::Open(m_path.c_str(), opt));
    if (!f || f->IsZombie())
        throw std::runtime_error("EventIO: failed to open " + m_path);
    f->Close();
}

std::string EventIO::sample_tree_name(const std::string &sample_name,
                                      const std::string &tree_prefix) const
{
    const std::string p = tree_prefix.empty() ? "events" : tree_prefix;
    return SnapshotService::sanitise_root_key(p) + "_" + SnapshotService::sanitise_root_key(sample_name);
}

ULong64_t EventIO::snapshot_event_list_merged(ROOT::RDF::RNode node,
                                              int sample_id,
                                              const std::string &sample_name,
                                              const std::vector<std::string> &columns,
                                              const std::string &selection,
                                              const std::string &tree_name_in) const
{
    return SnapshotService::snapshot_event_list_merged(std::move(node),
                                                       m_path,
                                                       sample_id,
                                                       sample_name,
                                                       columns,
                                                       selection,
                                                       tree_name_in);
}

ULong64_t EventIO::snapshot_event_list(ROOT::RDF::RNode node,
                                       const std::string &sample_name,
                                       const std::vector<std::string> &columns,
                                       const std::string &selection,
                                       const std::string &tree_prefix,
                                       bool overwrite_if_exists) const
{
    return SnapshotService::snapshot_event_list(std::move(node),
                                                m_path,
                                                sample_name,
                                                columns,
                                                selection,
                                                tree_prefix,
                                                overwrite_if_exists);
}

}
