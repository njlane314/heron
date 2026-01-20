/* -- C++ -- */
/**
 *  @file  io/src/AnalysisIO.cc
 *
 *  @brief Implementation for AnalysisIO helpers.
 */

#include "AnalysisIO.hh"

#include <ctime>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <TDirectory.h>
#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

namespace nuxsec
{
namespace io
{
namespace
{

std::string now_utc_iso8601()
{
    std::time_t t = std::time(nullptr);
    std::tm gm{};
#if defined(_WIN32)
    gmtime_s(&gm, &t);
#else
    gm = *std::gmtime(&t);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &gm);
    return std::string(buf);
}

const char *mode_to_string(AnalysisIO::OpenMode m)
{
    switch (m)
    {
    case AnalysisIO::OpenMode::kRead:
        return "READ";
    case AnalysisIO::OpenMode::kUpdate:
        return "UPDATE";
    default:
        return "READ";
    }
}

TDirectory *get_or_make_dir(TDirectory *parent, const std::string &name)
{
    if (!parent)
        return nullptr;
    if (auto *d = parent->GetDirectory(name.c_str()))
        return d;
    return parent->mkdir(name.c_str());
}

TDirectory *get_or_make_path(TDirectory *root, const std::string &path)
{
    TDirectory *cur = root;
    std::string token;
    std::istringstream is(path);
    while (std::getline(is, token, '/'))
    {
        if (token.empty())
            continue;
        cur = get_or_make_dir(cur, token);
        if (!cur)
            throw std::runtime_error("Failed to create directory path: " + path);
    }
    return cur;
}

TDirectory *get_dir_or_throw(TFile &f, const std::string &path)
{
    auto *d = f.GetDirectory(path.c_str());
    if (!d)
        throw std::runtime_error("Missing directory in analysis file: " + path);
    return d;
}

void write_objstring(TDirectory *dir, const std::string &key, const std::string &value)
{
    if (!dir)
        throw std::runtime_error("Null directory while writing meta key: " + key);
    dir->cd();
    TObjString s(value.c_str());
    s.Write(key.c_str(), TObject::kOverwrite);
}

std::string read_objstring(TDirectory *dir, const std::string &key)
{
    if (!dir)
        throw std::runtime_error("Null directory while reading meta key: " + key);
    auto *obj = dir->Get(key.c_str());
    auto *s = dynamic_cast<TObjString *>(obj);
    if (!s)
        throw std::runtime_error("Missing or wrong-type meta key: " + key);
    return std::string(s->String().Data());
}

} // namespace

class AnalysisIO::Impl
{
  public:
    std::string path;
    OpenMode mode;
    std::unique_ptr<TFile> file;

    mutable bool header_loaded = false;
    mutable AnalysisHeader header_cache;

    Impl(std::string p, OpenMode m) : path(std::move(p)), mode(m)
    {
        file = std::make_unique<TFile>(path.c_str(), mode_to_string(mode));
        if (!file || file->IsZombie())
            throw std::runtime_error("Failed to open analysis ROOT file: " + path);
    }

    void require_update() const
    {
        if (mode != OpenMode::kUpdate)
            throw std::runtime_error("AnalysisIO opened READ-only; reopen with UPDATE to write products");
    }

    const AnalysisHeader &load_header() const
    {
        if (header_loaded)
            return header_cache;

        auto *g = get_dir_or_throw(*file, "__global__");

        header_cache.schema = read_objstring(g, "schema");
        header_cache.analysis_name = read_objstring(g, "analysis_name");
        header_cache.analysis_tree = read_objstring(g, "analysis_tree");
        header_cache.created_utc = read_objstring(g, "created_utc");
        header_cache.sample_list_source = read_objstring(g, "sample_list_source");

        header_loaded = true;
        return header_cache;
    }
};

AnalysisIO::AnalysisIO(std::string analysis_root, OpenMode mode)
    : m(std::make_unique<Impl>(std::move(analysis_root), mode))
{
}

AnalysisIO::~AnalysisIO() = default;

AnalysisIO::AnalysisIO(AnalysisIO &&) noexcept = default;
AnalysisIO &AnalysisIO::operator=(AnalysisIO &&) noexcept = default;

const std::string &AnalysisIO::path() const noexcept { return m->path; }

const AnalysisHeader &AnalysisIO::header() const { return m->load_header(); }

std::string AnalysisIO::template_specs_1d_tsv() const
{
    auto *w = get_dir_or_throw(*m->file, "workspace");
    return read_objstring(w, "template_specs_1d_tsv");
}

std::vector<AnalysisSampleRef> AnalysisIO::samples() const
{
    auto *w = get_dir_or_throw(*m->file, "workspace");
    auto *t = dynamic_cast<TTree *>(w->Get("samples"));
    if (!t)
        throw std::runtime_error("Missing TTree workspace/samples in analysis file");

    std::vector<AnalysisSampleRef> out;

    std::string sample_name;
    std::string sample_rootio_path;
    int sample_kind = -1;
    int beam_mode = -1;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;

    t->SetBranchAddress("sample_name", &sample_name);
    t->SetBranchAddress("sample_rootio_path", &sample_rootio_path);
    t->SetBranchAddress("sample_kind", &sample_kind);
    t->SetBranchAddress("beam_mode", &beam_mode);
    t->SetBranchAddress("subrun_pot_sum", &subrun_pot_sum);
    t->SetBranchAddress("db_tortgt_pot_sum", &db_tortgt_pot_sum);
    t->SetBranchAddress("db_tor101_pot_sum", &db_tor101_pot_sum);

    const auto n = t->GetEntries();
    out.reserve(static_cast<size_t>(n));

    for (Long64_t i = 0; i < n; ++i)
    {
        t->GetEntry(i);
        AnalysisSampleRef s;
        s.sample_name = sample_name;
        s.sample_rootio_path = sample_rootio_path;
        s.sample_kind = sample_kind;
        s.beam_mode = beam_mode;
        s.subrun_pot_sum = subrun_pot_sum;
        s.db_tortgt_pot_sum = db_tortgt_pot_sum;
        s.db_tor101_pot_sum = db_tor101_pot_sum;
        out.push_back(std::move(s));
    }

    return out;
}

std::unique_ptr<TH1> AnalysisIO::get_hist1d(const std::string &family,
                                            const std::string &sample_name,
                                            const std::string &hist_name) const
{
    if (family.empty() || sample_name.empty() || hist_name.empty())
        throw std::runtime_error("AnalysisIO::get_hist1d: empty selector");

    const std::string dir_path = std::string("products/") + family + "/" + sample_name;
    auto *d = m->file->GetDirectory(dir_path.c_str());
    if (!d)
        throw std::runtime_error("Missing directory in analysis file: " + dir_path);

    auto *obj = d->Get(hist_name.c_str());
    auto *h = dynamic_cast<TH1 *>(obj);
    if (!h)
        throw std::runtime_error("Missing histogram: " + dir_path + "/" + hist_name);

    std::unique_ptr<TH1> out(dynamic_cast<TH1 *>(h->Clone(hist_name.c_str())));
    if (!out)
        throw std::runtime_error("Failed to clone histogram: " + hist_name);

    out->SetDirectory(nullptr);
    return out;
}

void AnalysisIO::put_histograms(const std::string &family,
                                const std::string &sample_name,
                                const std::vector<std::pair<std::string, const TH1 *>> &hists)
{
    m->require_update();

    if (family.empty())
        throw std::runtime_error("AnalysisIO::put_histograms: empty family");
    if (sample_name.empty())
        throw std::runtime_error("AnalysisIO::put_histograms: empty sample_name");

    const std::string dir_path = std::string("products/") + family + "/" + sample_name;
    TDirectory *d = get_or_make_path(m->file.get(), dir_path);
    d->cd();

    for (const auto &kv : hists)
    {
        const std::string &name = kv.first;
        const TH1 *h = kv.second;
        if (!h)
            continue;

        std::unique_ptr<TH1> c(dynamic_cast<TH1 *>(h->Clone(name.c_str())));
        if (!c)
            throw std::runtime_error("Failed to clone histogram: " + name);

        c->SetDirectory(d);
        c->Write(name.c_str(), TObject::kOverwrite);
    }
}

void AnalysisIO::flush()
{
    if (!m || !m->file)
        return;
    m->file->Write();
    m->file->Flush();
}

void AnalysisIO::init(const std::string &analysis_root,
                      const AnalysisHeader &header_in,
                      const std::vector<AnalysisSampleRef> &samples,
                      const std::string &template_specs_1d_tsv,
                      const std::string &template_specs_source)
{
    if (analysis_root.empty())
        throw std::runtime_error("AnalysisIO::init: empty output path");

    AnalysisHeader header = header_in;
    if (header.schema.empty())
        header.schema = kSchema;
    if (header.created_utc.empty())
        header.created_utc = now_utc_iso8601();

    TFile f(analysis_root.c_str(), "RECREATE");
    if (f.IsZombie())
        throw std::runtime_error("Failed to create analysis ROOT file: " + analysis_root);

    // __global__
    TDirectory *g = get_or_make_dir(&f, "__global__");
    write_objstring(g, "schema", header.schema);
    write_objstring(g, "analysis_name", header.analysis_name);
    write_objstring(g, "analysis_tree", header.analysis_tree);
    write_objstring(g, "created_utc", header.created_utc);
    write_objstring(g, "sample_list_source", header.sample_list_source);

    // workspace snapshots
    TDirectory *w = get_or_make_dir(&f, "workspace");
    write_objstring(w, "template_specs_1d_source", template_specs_source);
    write_objstring(w, "template_specs_1d_tsv", template_specs_1d_tsv);

    // workspace/samples index tree
    w->cd();
    TTree t("samples", "Nuxsec analysis sample index");

    std::string sample_name;
    std::string sample_rootio_path;
    int sample_kind = -1;
    int beam_mode = -1;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;

    t.Branch("sample_name", &sample_name);
    t.Branch("sample_rootio_path", &sample_rootio_path);
    t.Branch("sample_kind", &sample_kind);
    t.Branch("beam_mode", &beam_mode);
    t.Branch("subrun_pot_sum", &subrun_pot_sum);
    t.Branch("db_tortgt_pot_sum", &db_tortgt_pot_sum);
    t.Branch("db_tor101_pot_sum", &db_tor101_pot_sum);

    for (const auto &s : samples)
    {
        sample_name = s.sample_name;
        sample_rootio_path = s.sample_rootio_path;
        sample_kind = s.sample_kind;
        beam_mode = s.beam_mode;
        subrun_pot_sum = s.subrun_pot_sum;
        db_tortgt_pot_sum = s.db_tortgt_pot_sum;
        db_tor101_pot_sum = s.db_tor101_pot_sum;
        t.Fill();
    }

    t.Write("", TObject::kOverwrite);
    f.Write();
    f.Close();
}

} // namespace io
} // namespace nuxsec
