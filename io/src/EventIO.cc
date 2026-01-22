/* -- C++ -- */
/**
 *  @file  io/src/EventIO.cc
 *
 *  @brief Implementation of event-level IO.
 */

#include "EventIO.hh"

#include <cctype>
#include <memory>
#include <stdexcept>

#include <TArrayD.h>
#include <TArrayI.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TTree.h>

namespace nuxsec
{
namespace io
{

std::string EventIO::sanitise_root_key(std::string s)
{
    for (char &c : s)
    {
        const unsigned char u = static_cast<unsigned char>(c);
        if (!(std::isalnum(u) || c == '_'))
            c = '_';
    }
    if (s.empty())
        s = "sample";
    return s;
}

void EventIO::init(const std::string &out_path,
                   const EventHeader &header,
                   const std::vector<EventSampleRef> &sample_refs,
                   const std::string &event_schema_tsv,
                   const std::string &schema_tag)
{
    std::unique_ptr<TFile> fout(TFile::Open(out_path.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie())
        throw std::runtime_error("EventIO::init: failed to create output file: " + out_path);

    TObjString(header.analysis_name.c_str()).Write("analysis_name");
    TObjString(header.analysis_tree.c_str()).Write("analysis_tree");
    TObjString(header.sample_list_source.c_str()).Write("sample_list_source");

    if (!event_schema_tsv.empty())
    {
        const std::string key = schema_tag.empty() ? "event_schema"
                                                   : ("event_schema_" + sanitise_root_key(schema_tag));
        TObjString(event_schema_tsv.c_str()).Write(key.c_str());
    }

    TTree tref("sample_refs", "Sample references (source + POT/triggers)");
    std::string sample_name;
    std::string sample_rootio_path;
    int sample_kind = -1;
    int beam_mode = -1;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;

    tref.Branch("sample_name", &sample_name);
    tref.Branch("sample_rootio_path", &sample_rootio_path);
    tref.Branch("sample_kind", &sample_kind);
    tref.Branch("beam_mode", &beam_mode);
    tref.Branch("subrun_pot_sum", &subrun_pot_sum);
    tref.Branch("db_tortgt_pot_sum", &db_tortgt_pot_sum);
    tref.Branch("db_tor101_pot_sum", &db_tor101_pot_sum);

    for (const auto &r : sample_refs)
    {
        sample_name = r.sample_name;
        sample_rootio_path = r.sample_rootio_path;
        sample_kind = r.sample_kind;
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
    return sanitise_root_key(p) + "_" + sanitise_root_key(sample_name);
}

ULong64_t EventIO::snapshot_event_list(ROOT::RDF::RNode node,
                                       const std::string &sample_name,
                                       const std::string &x_column,
                                       const std::string &y_column,
                                       const std::vector<std::string> &double_columns,
                                       const std::vector<std::string> &int_columns,
                                       const std::string &selection,
                                       const std::string &tree_prefix,
                                       bool overwrite_if_exists) const
{
    ROOT::RDF::RNode filtered = std::move(node);
    if (!selection.empty() && selection != "true")
    {
        filtered = filtered.Filter(selection, "eventio_selection");
    }

    const std::string tree_name = sample_tree_name(sample_name, tree_prefix);

    std::vector<std::vector<double>> dvalues;
    dvalues.reserve(double_columns.size() + 2);
    dvalues.push_back(filtered.Take<double>(x_column).GetValue());
    if (!y_column.empty())
        dvalues.push_back(filtered.Take<double>(y_column).GetValue());
    for (const auto &col : double_columns)
        dvalues.push_back(filtered.Take<double>(col).GetValue());

    std::vector<std::vector<int>> ivalues;
    ivalues.reserve(int_columns.size());
    for (const auto &col : int_columns)
        ivalues.push_back(filtered.Take<int>(col).GetValue());

    const int dimensionality = y_column.empty() ? 1 : 2;
    const int ndoubles = static_cast<int>(dvalues.size());
    const int nints = static_cast<int>(ivalues.size());
    const int nevents = ndoubles > 0 ? static_cast<int>(dvalues.front().size()) : 0;

    for (const auto &vals : dvalues)
        if (static_cast<int>(vals.size()) != nevents)
            throw std::runtime_error("EventIO::snapshot_event_list: inconsistent double column sizes");
    for (const auto &vals : ivalues)
        if (static_cast<int>(vals.size()) != nevents)
            throw std::runtime_error("EventIO::snapshot_event_list: inconsistent int column sizes");

    TArrayD doubles(ndoubles * nevents);
    for (int v = 0; v < ndoubles; ++v)
        for (int i = 0; i < nevents; ++i)
            doubles[v * nevents + i] = dvalues[v][i];

    TArrayI integers(nints * nevents);
    for (int v = 0; v < nints; ++v)
        for (int i = 0; i < nevents; ++i)
            integers[v * nevents + i] = ivalues[v][i];

    TObjArray value_names;
    int name_index = 0;
    value_names.AddAtAndExpand(new TObjString("X"), name_index++);
    if (dimensionality > 1)
        value_names.AddAtAndExpand(new TObjString("Y"), name_index++);
    for (const auto &name : double_columns)
        value_names.AddAtAndExpand(new TObjString(name.c_str()), name_index++);
    for (const auto &name : int_columns)
        value_names.AddAtAndExpand(new TObjString(name.c_str()), name_index++);

    std::unique_ptr<TFile> fout(TFile::Open(m_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
        throw std::runtime_error("EventIO::snapshot_event_list: failed to open " + m_path);

    if (overwrite_if_exists)
        fout->Delete((tree_name + ";*").c_str());

    TTree tree(tree_name.c_str(), "Event list");
    int dim_i = dimensionality;
    int ndoubles_i = ndoubles;
    int nints_i = nints;
    int nevents_i = nevents;
    tree.Branch("dimensionality", &dim_i, "dimensionality/I");
    tree.Branch("ndoubles", &ndoubles_i, "ndoubles/I");
    tree.Branch("nints", &nints_i, "nints/I");
    tree.Branch("nevents", &nevents_i, "nevents/I");
    tree.Branch("doubles", &doubles);
    tree.Branch("integers", &integers);
    tree.Branch("value_names", &value_names);
    tree.Fill();
    tree.Write(tree_name.c_str(), TObject::kOverwrite);

    fout->Close();
    return static_cast<ULong64_t>(nevents);
}

} // namespace io
} // namespace nuxsec
