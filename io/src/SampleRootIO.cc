/* -- C++ -- */
/**
 *  @file  io/src/SampleRootIO.cc
 *
 *  @brief Implementation for SampleRootIO helpers.
 */

#include "SampleRootIO.hh"

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <memory>
#include <stdexcept>
#include <utility>

namespace nuxsec
{

void SampleRootIO::write(const Sample &sample, const std::string &out_file)
{
    std::unique_ptr<TFile> f(TFile::Open(out_file.c_str(), "UPDATE"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open merged output file for UPDATE: " + out_file);
    }

    TDirectory *d = f->GetDirectory("SampleRootIO");
    if (!d)
    {
        d = f->mkdir("SampleRootIO");
    }
    d->cd();

    TNamed("sample_name", sample.sample_name.c_str()).Write("sample_name", TObject::kOverwrite);
    TNamed("sample_kind", sample_kind_name(sample.kind)).Write("sample_kind", TObject::kOverwrite);
    TNamed("beam_mode", beam_mode_name(sample.beam)).Write("beam_mode", TObject::kOverwrite);

    TParameter<double>("subrun_pot_sum", sample.subrun_pot_sum).Write("subrun_pot_sum", TObject::kOverwrite);
    TParameter<double>("db_tortgt_pot_sum", sample.db_tortgt_pot_sum)
        .Write("db_tortgt_pot_sum", TObject::kOverwrite);
    TParameter<double>("db_tor101_pot_sum", sample.db_tor101_pot_sum)
        .Write("db_tor101_pot_sum", TObject::kOverwrite);
    TParameter<double>("normalisation", sample.normalisation).Write("normalisation", TObject::kOverwrite);
    TParameter<double>("normalised_pot_sum", sample.normalised_pot_sum)
        .Write("normalised_pot_sum", TObject::kOverwrite);

    {
        TTree fragments("fragments", "Art file fragments included in sample aggregation");

        std::string fragment_name;
        std::string artio_path;
        double subrun_pot_sum = 0.0;
        double db_tortgt_pot = 0.0;
        double db_tor101_pot = 0.0;
        double normalisation = 1.0;
        double normalised_pot_sum = 0.0;

        fragments.Branch("fragment_name", &fragment_name);
        fragments.Branch("artio_path", &artio_path);
        fragments.Branch("subrun_pot_sum", &subrun_pot_sum);
        fragments.Branch("db_tortgt_pot", &db_tortgt_pot);
        fragments.Branch("db_tor101_pot", &db_tor101_pot);
        fragments.Branch("normalisation", &normalisation);
        fragments.Branch("normalised_pot_sum", &normalised_pot_sum);

        for (const auto &fragment : sample.fragments)
        {
            fragment_name = fragment.fragment_name;
            artio_path = fragment.artio_path;
            subrun_pot_sum = fragment.subrun_pot_sum;
            db_tortgt_pot = fragment.db_tortgt_pot;
            db_tor101_pot = fragment.db_tor101_pot;
            normalisation = fragment.normalisation;
            normalised_pot_sum = fragment.normalised_pot_sum;
            fragments.Fill();
        }

        fragments.Write("fragments", TObject::kOverwrite);
    }

    f->Write();
    f->Close();
}

Sample SampleRootIO::read(const std::string &in_file)
{
    std::unique_ptr<TFile> f(TFile::Open(in_file.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open merged output file for READ: " + in_file);
    }

    TDirectory *d = f->GetDirectory("SampleRootIO");
    if (!d)
    {
        throw std::runtime_error("Missing SampleRootIO directory in file: " + in_file);
    }
    d->cd();

    Sample out;
    {
        TObject *obj = d->Get("sample_name");
        auto *named = dynamic_cast<TNamed *>(obj);
        if (!named)
        {
            throw std::runtime_error("Missing sample_name metadata in SampleRootIO directory");
        }
        out.sample_name = named->GetTitle();
    }

    {
        TObject *obj = d->Get("sample_kind");
        auto *named = dynamic_cast<TNamed *>(obj);
        if (!named)
        {
            throw std::runtime_error("Missing sample_kind metadata in SampleRootIO directory");
        }
        out.kind = parse_sample_kind(named->GetTitle());
    }

    {
        TObject *obj = d->Get("beam_mode");
        auto *named = dynamic_cast<TNamed *>(obj);
        if (!named)
        {
            throw std::runtime_error("Missing beam_mode metadata in SampleRootIO directory");
        }
        out.beam = parse_beam_mode(named->GetTitle());
    }

    auto read_param_double = [d](const char *key)
    {
        TObject *obj = d->Get(key);
        auto *param = dynamic_cast<TParameter<double> *>(obj);
        if (!param)
        {
            throw std::runtime_error("Missing TParameter<double> for key: " + std::string(key));
        }
        return param->GetVal();
    };

    out.subrun_pot_sum = read_param_double("subrun_pot_sum");
    out.db_tortgt_pot_sum = read_param_double("db_tortgt_pot_sum");
    out.db_tor101_pot_sum = read_param_double("db_tor101_pot_sum");
    out.normalisation = read_param_double("normalisation");
    out.normalised_pot_sum = read_param_double("normalised_pot_sum");

    TObject *obj = d->Get("fragments");
    auto *tree = dynamic_cast<TTree *>(obj);
    if (!tree)
    {
        throw std::runtime_error("Missing fragments tree in SampleRootIO directory");
    }

    std::string fragment_name;
    std::string artio_path;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot = 0.0;
    double db_tor101_pot = 0.0;
    double normalisation = 1.0;
    double normalised_pot_sum = 0.0;

    tree->SetBranchAddress("fragment_name", &fragment_name);
    tree->SetBranchAddress("artio_path", &artio_path);
    tree->SetBranchAddress("subrun_pot_sum", &subrun_pot_sum);
    tree->SetBranchAddress("db_tortgt_pot", &db_tortgt_pot);
    tree->SetBranchAddress("db_tor101_pot", &db_tor101_pot);
    tree->SetBranchAddress("normalisation", &normalisation);
    tree->SetBranchAddress("normalised_pot_sum", &normalised_pot_sum);

    const Long64_t n = tree->GetEntries();
    out.fragments.reserve(static_cast<size_t>(n));
    for (Long64_t i = 0; i < n; ++i)
    {
        tree->GetEntry(i);
        SampleFragment fragment;
        fragment.fragment_name = fragment_name;
        fragment.artio_path = artio_path;
        fragment.subrun_pot_sum = subrun_pot_sum;
        fragment.db_tortgt_pot = db_tortgt_pot;
        fragment.db_tor101_pot = db_tor101_pot;
        fragment.normalisation = normalisation;
        fragment.normalised_pot_sum = normalised_pot_sum;
        out.fragments.push_back(std::move(fragment));
    }

    return out;
}

} // namespace nuxsec
