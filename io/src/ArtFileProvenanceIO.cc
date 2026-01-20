/* -- C++ -- */
/**
 *  @file  io/src/ArtFileProvenanceIO.cc
 *
 *  @brief Implementation for Art file provenance ROOT IO.
 */

#include "ArtFileProvenanceIO.hh"

namespace nuxsec
{

void ArtFileProvenanceIO::write(const artio::Provenance &r, const std::string &out_file)
{
    std::unique_ptr<TFile> f(TFile::Open(out_file.c_str(), "UPDATE"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open merged output file for UPDATE: " + out_file);
    }

    TDirectory *d = f->GetDirectory("nuxsec_art_provenance");
    if (!d)
    {
        d = f->mkdir("nuxsec_art_provenance");
    }
    d->cd();

    TNamed("stage_name", r.cfg.stage_name.c_str()).Write("stage_name", TObject::kOverwrite);
    TNamed("sample_kind", sample_kind_name(r.kind)).Write("sample_kind", TObject::kOverwrite);
    TNamed("beam_mode", beam_mode_name(r.beam)).Write("beam_mode", TObject::kOverwrite);

    TParameter<double>("subrun_pot_sum", r.subrun.pot_sum).Write("subrun_pot_sum", TObject::kOverwrite);
    TParameter<long long>("subrun_entries", r.subrun.n_entries).Write("subrun_entries", TObject::kOverwrite);
    TParameter<long long>("unique_run_subrun_pairs", static_cast<long long>(r.subrun.unique_pairs.size()))
        .Write("unique_run_subrun_pairs", TObject::kOverwrite);

    TParameter<double>("scale_factor", r.scale).Write("scale_factor", TObject::kOverwrite);

    {
        TObjArray arr;
        arr.SetOwner(true);
        for (const auto &in : r.input_files)
        {
            arr.Add(new TObjString(in.c_str()));
        }
        arr.Write("input_files", TObject::kSingleKey | TObject::kOverwrite);
    }

    {
        TTree rs("run_subrun", "Unique (run,subrun) pairs used for DB sums");
        Int_t run = 0;
        Int_t subrun = 0;
        rs.Branch("run", &run, "run/I");
        rs.Branch("subrun", &subrun, "subrun/I");
        for (const auto &p : r.subrun.unique_pairs)
        {
            run = p.run;
            subrun = p.subrun;
            rs.Fill();
        }
        rs.Write("run_subrun", TObject::kOverwrite);
    }

    f->Write();
    f->Close();
}

artio::Provenance ArtFileProvenanceIO::read(const std::string &in_file)
{
    std::unique_ptr<TFile> f(TFile::Open(in_file.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open merged output file for READ: " + in_file);
    }

    TDirectory *d = f->GetDirectory("nuxsec_art_provenance");
    if (!d)
    {
        throw std::runtime_error("Missing nuxsec_art_provenance directory in file: " + in_file);
    }
    d->cd();

    const SampleKind kind = parse_sample_kind(read_named_string(d, "sample_kind"));
    const BeamMode beam = parse_beam_mode(read_named_string(d, "beam_mode"));

    return read_directory(d, kind, beam);
}

artio::Provenance ArtFileProvenanceIO::read(const std::string &in_file, SampleKind kind, BeamMode beam)
{
    std::unique_ptr<TFile> f(TFile::Open(in_file.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open merged output file for READ: " + in_file);
    }

    TDirectory *d = f->GetDirectory("nuxsec_art_provenance");
    if (!d)
    {
        throw std::runtime_error("Missing nuxsec_art_provenance directory in file: " + in_file);
    }
    d->cd();

    return read_directory(d, kind, beam);
}

artio::Provenance ArtFileProvenanceIO::read_directory(TDirectory *d, SampleKind kind, BeamMode beam)
{
    artio::Provenance r;
    r.cfg.stage_name = read_named_string(d, "stage_name");
    r.kind = kind;
    r.beam = beam;

    r.subrun.pot_sum = read_param<double>(d, "subrun_pot_sum");
    r.subrun.n_entries = read_param<long long>(d, "subrun_entries");

    r.scale = read_param<double>(d, "scale_factor");

    r.input_files = read_input_files(d);
    r.subrun.unique_pairs = read_run_subrun_pairs(d);
    return r;
}

std::string ArtFileProvenanceIO::read_named_string(TDirectory *d, const char *key)
{
    TObject *obj = d->Get(key);
    auto *named = dynamic_cast<TNamed *>(obj);
    if (!named)
    {
        throw std::runtime_error("Missing TNamed for key: " + std::string(key));
    }
    return std::string(named->GetTitle());
}

std::vector<std::string> ArtFileProvenanceIO::read_input_files(TDirectory *d)
{
    std::vector<std::string> files;
    TObject *obj = d->Get("input_files");
    auto *arr = dynamic_cast<TObjArray *>(obj);
    if (!arr)
    {
        throw std::runtime_error("Missing input_files array");
    }
    const int n = arr->GetEntries();
    files.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        auto *entry = dynamic_cast<TObjString *>(arr->At(i));
        if (!entry)
        {
            throw std::runtime_error("Invalid entry in input_files array");
        }
        files.emplace_back(entry->GetString().Data());
    }
    return files;
}

std::vector<artio::RunSubrunPair> ArtFileProvenanceIO::read_run_subrun_pairs(TDirectory *d)
{
    TObject *obj = d->Get("run_subrun");
    auto *tree = dynamic_cast<TTree *>(obj);
    if (!tree)
    {
        throw std::runtime_error("Missing run_subrun tree");
    }

    Int_t run = 0;
    Int_t subrun = 0;
    tree->SetBranchAddress("run", &run);
    tree->SetBranchAddress("subrun", &subrun);

    const Long64_t n = tree->GetEntries();
    std::vector<artio::RunSubrunPair> pairs;
    pairs.reserve(static_cast<size_t>(n));
    for (Long64_t i = 0; i < n; ++i)
    {
        tree->GetEntry(i);
        pairs.push_back(artio::RunSubrunPair{static_cast<int>(run), static_cast<int>(subrun)});
    }
    return pairs;
}

} // namespace nuxsec
