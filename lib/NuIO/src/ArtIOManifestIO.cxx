/**
 *  @file  lib/NuIO/src/ArtIOManifestIO.cxx
 *
 *  @brief Implementation of ArtIO manifest metadata IO
 */

#include "NuIO/ArtIOManifestIO.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TObject.h>
#include <TParameter.h>
#include <TTree.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace nuio
{

std::vector<std::string> ArtIOManifestIO::ListStages(const std::string &artio_file)
{
    std::unique_ptr<TFile> f(TFile::Open(artio_file.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
        return {};
    }

    auto *t = dynamic_cast<TTree *>(f->Get("Stages"));
    if (!t)
    {
        return {};
    }

    std::string stage_name;
    t->SetBranchAddress("stage_name", &stage_name);

    const Long64_t n = t->GetEntries();
    std::vector<std::string> out;
    out.reserve(static_cast<size_t>(n));
    for (Long64_t i = 0; i < n; ++i)
    {
        t->GetEntry(i);
        out.push_back(stage_name);
    }
    return out;
}

static TDirectory *GetOrMakeDir(TFile &f, const char *name)
{
    TDirectory *d = f.GetDirectory(name);
    if (!d)
    {
        d = f.mkdir(name);
    }
    if (!d)
    {
        throw std::runtime_error(std::string("Failed to create directory: ") + name);
    }
    return d;
}

void ArtIOManifestIO::WriteStage(const std::string &artio_file,
                                 const std::string &db_path,
                                 double pot_scale,
                                 const ArtProvenance &stage)
{
    std::unique_ptr<TFile> f(TFile::Open(artio_file.c_str(), "UPDATE"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open ArtIO file for UPDATE: " + artio_file);
    }

    std::string stage_name;
    std::string filelist_path;
    std::string kind;
    std::string beam;
    Long64_t n_input_files = 0;
    Double_t subrun_pot_sum = 0.0;
    Long64_t subrun_entries = 0;
    Long64_t n_unique_pairs = 0;

    Double_t tortgt_sum = 0.0;
    Double_t tor101_sum = 0.0;
    Double_t tor860_sum = 0.0;
    Double_t tor875_sum = 0.0;
    Long64_t EA9CNT_sum = 0;
    Long64_t E1DCNT_sum = 0;
    Long64_t EXTTrig_sum = 0;
    Long64_t Gate1Trig_sum = 0;
    Long64_t Gate2Trig_sum = 0;

    auto *tStages = new TTree("Stages", "ArtIO stage inventory");
    tStages->Branch("stage_name", &stage_name);
    tStages->Branch("filelist_path", &filelist_path);
    tStages->Branch("kind", &kind);
    tStages->Branch("beam", &beam);
    tStages->Branch("n_input_files", &n_input_files);
    tStages->Branch("subrun_pot_sum", &subrun_pot_sum);
    tStages->Branch("subrun_entries", &subrun_entries);
    tStages->Branch("n_unique_pairs", &n_unique_pairs);
    tStages->Branch("tortgt_sum", &tortgt_sum);
    tStages->Branch("tor101_sum", &tor101_sum);
    tStages->Branch("tor860_sum", &tor860_sum);
    tStages->Branch("tor875_sum", &tor875_sum);
    tStages->Branch("EA9CNT_sum", &EA9CNT_sum);
    tStages->Branch("E1DCNT_sum", &E1DCNT_sum);
    tStages->Branch("EXTTrig_sum", &EXTTrig_sum);
    tStages->Branch("Gate1Trig_sum", &Gate1Trig_sum);
    tStages->Branch("Gate2Trig_sum", &Gate2Trig_sum);

    std::string p_stage_name;
    Int_t p_run = 0;
    Int_t p_subrun = 0;

    auto *tPairs = new TTree("RunSubruns", "Run/subrun inventory keyed by stage_name");
    tPairs->Branch("stage_name", &p_stage_name);
    tPairs->Branch("run", &p_run, "run/I");
    tPairs->Branch("subrun", &p_subrun, "subrun/I");

    stage_name = stage.cfg.stage_name;
    filelist_path = stage.cfg.filelist_path;
    kind = SampleKindName(stage.kind);
    beam = BeamModeName(stage.beam);
    n_input_files = static_cast<Long64_t>(stage.input_files.size());

    subrun_pot_sum = stage.subrun.pot_sum;
    subrun_entries = static_cast<Long64_t>(stage.subrun.n_entries);
    n_unique_pairs = static_cast<Long64_t>(stage.subrun.unique_pairs.size());

    tortgt_sum = stage.runinfo.tortgt_sum;
    tor101_sum = stage.runinfo.tor101_sum;
    tor860_sum = stage.runinfo.tor860_sum;
    tor875_sum = stage.runinfo.tor875_sum;
    EA9CNT_sum = static_cast<Long64_t>(stage.runinfo.EA9CNT_sum);
    E1DCNT_sum = static_cast<Long64_t>(stage.runinfo.E1DCNT_sum);
    EXTTrig_sum = static_cast<Long64_t>(stage.runinfo.EXTTrig_sum);
    Gate1Trig_sum = static_cast<Long64_t>(stage.runinfo.Gate1Trig_sum);
    Gate2Trig_sum = static_cast<Long64_t>(stage.runinfo.Gate2Trig_sum);

    tStages->Fill();

    p_stage_name = stage.cfg.stage_name;
    for (const auto &rs : stage.subrun.unique_pairs)
    {
        p_run = rs.run;
        p_subrun = rs.subrun;
        tPairs->Fill();
    }

    f->cd();
    tStages->Write("Stages", TObject::kOverwrite);
    tPairs->Write("RunSubruns", TObject::kOverwrite);

    TDirectory *d = GetOrMakeDir(*f, "ArtIO");
    d->cd();
    TNamed("db_path", db_path.c_str()).Write("db_path", TObject::kOverwrite);
    TParameter<double>("pot_scale", pot_scale).Write("pot_scale", TObject::kOverwrite);

    f->Write();
    f->Close();
}

}
