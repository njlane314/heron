#include "EventListIO.hh"

#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

#include "PlottingHelper.hh"
#include "SampleIO.hh"
#include "SnapshotService.hh"

namespace
{
std::string read_objstring_optional(TFile &f, const char *key)
{
    TObject *obj = f.Get(key);
    auto *s = dynamic_cast<TObjString *>(obj);
    if (!s)
        return {};
    return s->GetString().Data();
}

}

namespace nu
{
void EventListIO::init(const std::string &out_path,
                       const EventListHeader &header,
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
                "EventListIO::init: failed to create output directory: "
                + output_path.parent_path().string()
                + " (" + ec.message() + ")");
        }
    }

    std::unique_ptr<TFile> fout(TFile::Open(out_path.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie())
        throw std::runtime_error("EventListIO::init: failed to create output file: " + out_path);

    TObjString(header.analysis_name.c_str()).Write("analysis_name");
    TObjString(header.provenance_tree.c_str()).Write("provenance_tree");
    TObjString(header.event_tree.c_str()).Write("event_tree");
    TObjString(header.sample_list_source.c_str()).Write("sample_list_source");
    if (!header.heron_set.empty())
    {
        TObjString(header.heron_set.c_str()).Write("heron_set");
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

EventListIO EventListIO::read(std::string path)
{
    return EventListIO(std::move(path), OpenMode::kRead);
}

EventListIO::EventListIO(std::string path, OpenMode mode)
    : m_path(std::move(path)), m_mode(mode)
{
    const char *opt = (m_mode == OpenMode::kRead) ? "READ" : "UPDATE";
    std::unique_ptr<TFile> fin(TFile::Open(m_path.c_str(), opt));
    if (!fin || fin->IsZombie())
        throw std::runtime_error("EventListIO: failed to open " + m_path);

    m_header.analysis_name = read_objstring_optional(*fin, "analysis_name");
    m_header.provenance_tree = read_objstring_optional(*fin, "provenance_tree");
    m_header.event_tree = read_objstring_optional(*fin, "event_tree");
    m_header.sample_list_source = read_objstring_optional(*fin, "sample_list_source");
    m_header.heron_set = read_objstring_optional(*fin, "heron_set");
    m_header.event_output_dir = read_objstring_optional(*fin, "event_output_dir");

    auto *t = dynamic_cast<TTree *>(fin->Get("sample_refs"));
    if (!t)
        throw std::runtime_error("EventListIO: missing sample_refs tree in " + m_path);

    int sample_id = -1;
    std::string *sample_name = nullptr;
    std::string *sample_rootio_path = nullptr;
    int sample_origin = -1;
    int beam_mode = -1;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;

    t->SetBranchAddress("sample_id", &sample_id);
    t->SetBranchAddress("sample_name", &sample_name);
    t->SetBranchAddress("sample_rootio_path", &sample_rootio_path);
    t->SetBranchAddress("sample_origin", &sample_origin);
    t->SetBranchAddress("beam_mode", &beam_mode);
    t->SetBranchAddress("subrun_pot_sum", &subrun_pot_sum);
    t->SetBranchAddress("db_tortgt_pot_sum", &db_tortgt_pot_sum);
    t->SetBranchAddress("db_tor101_pot_sum", &db_tor101_pot_sum);

    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i)
    {
        t->GetEntry(i);
        if (!sample_name || !sample_rootio_path)
            throw std::runtime_error("EventListIO: sample_refs missing string branches");

        SampleInfo info;
        info.sample_name = *sample_name;
        info.sample_rootio_path = *sample_rootio_path;
        info.sample_origin = sample_origin;
        info.beam_mode = beam_mode;
        info.subrun_pot_sum = subrun_pot_sum;
        info.db_tortgt_pot_sum = db_tortgt_pot_sum;
        info.db_tor101_pot_sum = db_tor101_pot_sum;

        m_sample_refs.emplace(sample_id, std::move(info));
        if (sample_id > m_max_sample_id)
            m_max_sample_id = sample_id;
    }

    fin->Close();
}

std::string EventListIO::event_tree() const
{
    return m_header.event_tree.empty() ? "events" : m_header.event_tree;
}

std::string EventListIO::sample_tree_name(const std::string &sample_name,
                                          const std::string &tree_prefix) const
{
    const std::string p = tree_prefix.empty() ? "events" : tree_prefix;
    return SnapshotService::sanitise_root_key(p) + "_" + SnapshotService::sanitise_root_key(sample_name);
}

ULong64_t EventListIO::snapshot_event_list_merged(ROOT::RDF::RNode node,
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

ULong64_t EventListIO::snapshot_event_list(ROOT::RDF::RNode node,
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

ROOT::RDataFrame EventListIO::rdf() const
{
    return ROOT::RDataFrame(event_tree(), m_path);
}

std::shared_ptr<const std::vector<char>> EventListIO::mask_for_origin(SampleIO::SampleOrigin origin) const
{
    const int want = static_cast<int>(origin);
    auto mask = std::make_shared<std::vector<char>>(static_cast<size_t>(m_max_sample_id + 1), 0);

    for (const auto &kv : m_sample_refs)
    {
        const int sid = kv.first;
        const SampleInfo &info = kv.second;
        if (sid >= 0 && sid <= m_max_sample_id)
            (*mask)[static_cast<size_t>(sid)] = (info.sample_origin == want) ? 1 : 0;
    }
    return mask;
}

std::shared_ptr<const std::vector<char>> EventListIO::mask_for_data() const
{
    return mask_for_origin(SampleIO::SampleOrigin::kData);
}

std::shared_ptr<const std::vector<char>> EventListIO::mask_for_ext() const
{
    return mask_for_origin(SampleIO::SampleOrigin::kEXT);
}

std::shared_ptr<const std::vector<char>> EventListIO::mask_for_mc_like() const
{
    auto mask = std::make_shared<std::vector<char>>(static_cast<size_t>(m_max_sample_id + 1), 0);
    const int data_id = static_cast<int>(SampleIO::SampleOrigin::kData);

    for (const auto &kv : m_sample_refs)
    {
        const int sid = kv.first;
        const SampleInfo &info = kv.second;
        if (sid < 0 || sid > m_max_sample_id)
            continue;
        (*mask)[static_cast<size_t>(sid)] = (info.sample_origin != data_id) ? 1 : 0;
    }
    return mask;
}

double EventListIO::total_pot_data() const
{
    const int data_id = static_cast<int>(SampleIO::SampleOrigin::kData);
    double tot = 0.0;
    for (const auto &kv : m_sample_refs)
    {
        const SampleInfo &info = kv.second;
        if (info.sample_origin == data_id)
            tot += info.db_tortgt_pot_sum;
    }
    return tot;
}

double EventListIO::total_pot_mc() const
{
    const int data_id = static_cast<int>(SampleIO::SampleOrigin::kData);
    double tot = 0.0;
    for (const auto &kv : m_sample_refs)
    {
        const SampleInfo &info = kv.second;
        if (info.sample_origin != data_id)
            tot += info.db_tortgt_pot_sum;
    }
    return tot;
}

std::string EventListIO::beamline_label() const
{
    int seen = -1;
    for (const auto &kv : m_sample_refs)
    {
        const int b = kv.second.beam_mode;
        if (b < 0)
            continue;
        if (seen < 0)
            seen = b;
        else if (b != seen)
            return "mixed";
    }

    if (seen == static_cast<int>(SampleIO::BeamMode::kNuMI))
        return "numi";
    if (seen == static_cast<int>(SampleIO::BeamMode::kBNB))
        return "bnb";
    return "unknown";
}

}
