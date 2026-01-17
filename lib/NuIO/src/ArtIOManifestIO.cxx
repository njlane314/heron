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
        return {};

    auto *t = dynamic_cast<TTree *>(f->Get("Stages"));
    if (!t)
        return {};

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
        d = f.mkdir(name);
    if (!d)
        throw std::runtime_error(std::string("Failed to create directory: ") + name);
    return d;
}

void ArtIOManifestIO::AppendStages(const std::string &artio_file,
                                   const std::string &db_path,
                                   double pot_scale,
                                   const std::vector<ArtIOStage> &stages)
{
    if (stages.empty())
        return;

    std::unique_ptr<TFile> f(TFile::Open(artio_file.c_str(), "UPDATE"));
    if (!f || f->IsZombie())
        throw std::runtime_error("Failed to open ArtIO file for UPDATE: " + artio_file);

    auto *tStages = dynamic_cast<TTree *>(f->Get("Stages"));
    auto *tPairs = dynamic_cast<TTree *>(f->Get("RunSubruns"));

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

    const bool newStages = (tStages == nullptr);
    if (newStages)
    {
        tStages = new TTree("Stages", "ArtIO stage inventory");
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
    }
    else
    {
        tStages->SetBranchAddress("stage_name", &stage_name);
        tStages->SetBranchAddress("filelist_path", &filelist_path);
        tStages->SetBranchAddress("kind", &kind);
        tStages->SetBranchAddress("beam", &beam);
        tStages->SetBranchAddress("n_input_files", &n_input_files);
        tStages->SetBranchAddress("subrun_pot_sum", &subrun_pot_sum);
        tStages->SetBranchAddress("subrun_entries", &subrun_entries);
        tStages->SetBranchAddress("n_unique_pairs", &n_unique_pairs);
        tStages->SetBranchAddress("tortgt_sum", &tortgt_sum);
        tStages->SetBranchAddress("tor101_sum", &tor101_sum);
        tStages->SetBranchAddress("tor860_sum", &tor860_sum);
        tStages->SetBranchAddress("tor875_sum", &tor875_sum);
        tStages->SetBranchAddress("EA9CNT_sum", &EA9CNT_sum);
        tStages->SetBranchAddress("E1DCNT_sum", &E1DCNT_sum);
        tStages->SetBranchAddress("EXTTrig_sum", &EXTTrig_sum);
        tStages->SetBranchAddress("Gate1Trig_sum", &Gate1Trig_sum);
        tStages->SetBranchAddress("Gate2Trig_sum", &Gate2Trig_sum);
    }

    std::string p_stage_name;
    Int_t p_run = 0;
    Int_t p_subrun = 0;

    const bool newPairs = (tPairs == nullptr);
    if (newPairs)
    {
        tPairs = new TTree("RunSubruns", "Run/subrun inventory keyed by stage_name");
        tPairs->Branch("stage_name", &p_stage_name);
        tPairs->Branch("run", &p_run, "run/I");
        tPairs->Branch("subrun", &p_subrun, "subrun/I");
    }
    else
    {
        tPairs->SetBranchAddress("stage_name", &p_stage_name);
        tPairs->SetBranchAddress("run", &p_run);
        tPairs->SetBranchAddress("subrun", &p_subrun);
    }

    for (const auto &s : stages)
    {
        stage_name = s.cfg.stage_name;
        filelist_path = s.cfg.filelist_path;
        kind = SampleKindName(s.kind);
        beam = BeamModeName(s.beam);
        n_input_files = static_cast<Long64_t>(s.n_input_files);

        subrun_pot_sum = s.subrun.pot_sum;
        subrun_entries = static_cast<Long64_t>(s.subrun.n_entries);
        n_unique_pairs = static_cast<Long64_t>(s.subrun.unique_pairs.size());

        tortgt_sum = s.runinfo.tortgt_sum;
        tor101_sum = s.runinfo.tor101_sum;
        tor860_sum = s.runinfo.tor860_sum;
        tor875_sum = s.runinfo.tor875_sum;
        EA9CNT_sum = static_cast<Long64_t>(s.runinfo.EA9CNT_sum);
        E1DCNT_sum = static_cast<Long64_t>(s.runinfo.E1DCNT_sum);
        EXTTrig_sum = static_cast<Long64_t>(s.runinfo.EXTTrig_sum);
        Gate1Trig_sum = static_cast<Long64_t>(s.runinfo.Gate1Trig_sum);
        Gate2Trig_sum = static_cast<Long64_t>(s.runinfo.Gate2Trig_sum);

        tStages->Fill();

        p_stage_name = s.cfg.stage_name;
        for (const auto &rs : s.subrun.unique_pairs)
        {
            p_run = rs.run;
            p_subrun = rs.subrun;
            tPairs->Fill();
        }
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
