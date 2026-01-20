/* -- C++ -- */
/**
 *  @file  io/src/SampleIO.cc
 *
 *  @brief Implementation for SampleIO helpers.
 */

#include "SampleIO.hh"

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <utility>

namespace nuxsec
{

namespace sample
{

const char *SampleIO::sample_kind_name(SampleKind k)
{
    switch (k)
    {
    case SampleKind::kData:
        return "data";
    case SampleKind::kEXT:
        return "ext";
    case SampleKind::kOverlay:
        return "mc_overlay";
    case SampleKind::kDirt:
        return "mc_dirt";
    case SampleKind::kStrangeness:
        return "mc_strangeness";
    default:
        return "unknown";
    }
}

SampleIO::SampleKind SampleIO::parse_sample_kind(const std::string &name)
{
    std::string lowered = name;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                   [](unsigned char c)
                   {
                       return static_cast<char>(std::tolower(c));
                   });

    if (lowered == "data")
    {
        return SampleKind::kData;
    }
    if (lowered == "ext")
    {
        return SampleKind::kEXT;
    }
    if (lowered == "overlay")
    {
        return SampleKind::kOverlay;
    }
    if (lowered == "dirt")
    {
        return SampleKind::kDirt;
    }
    if (lowered == "strangeness")
    {
        return SampleKind::kStrangeness;
    }
    return SampleKind::kUnknown;
}

const char *SampleIO::beam_mode_name(BeamMode b)
{
    switch (b)
    {
    case BeamMode::kNuMI:
        return "numi";
    case BeamMode::kBNB:
        return "bnb";
    default:
        return "unknown";
    }
}

SampleIO::BeamMode SampleIO::parse_beam_mode(const std::string &name)
{
    std::string lowered = name;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                   [](unsigned char c)
                   {
                       return static_cast<char>(std::tolower(c));
                   });

    if (lowered == "numi")
    {
        return BeamMode::kNuMI;
    }
    if (lowered == "bnb")
    {
        return BeamMode::kBNB;
    }
    return BeamMode::kUnknown;
}

void SampleIO::write(const Sample &sample, const std::string &out_file)
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
        TTree entries("entries", "Art file entries included in sample aggregation");

        std::string entry_name;
        std::string artio_path;
        double subrun_pot_sum = 0.0;
        double db_tortgt_pot = 0.0;
        double db_tor101_pot = 0.0;
        double normalisation = 1.0;
        double normalised_pot_sum = 0.0;

        entries.Branch("entry_name", &entry_name);
        entries.Branch("artio_path", &artio_path);
        entries.Branch("subrun_pot_sum", &subrun_pot_sum);
        entries.Branch("db_tortgt_pot", &db_tortgt_pot);
        entries.Branch("db_tor101_pot", &db_tor101_pot);
        entries.Branch("normalisation", &normalisation);
        entries.Branch("normalised_pot_sum", &normalised_pot_sum);

        for (const auto &entry : sample.entries)
        {
            entry_name = entry.entry_name;
            artio_path = entry.artio_path;
            subrun_pot_sum = entry.subrun_pot_sum;
            db_tortgt_pot = entry.db_tortgt_pot;
            db_tor101_pot = entry.db_tor101_pot;
            normalisation = entry.normalisation;
            normalised_pot_sum = entry.normalised_pot_sum;
            entries.Fill();
        }

        entries.Write("entries", TObject::kOverwrite);
    }

    f->Write();
    f->Close();
}

SampleIO::Sample SampleIO::read(const std::string &in_file)
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

    TObject *obj = d->Get("entries");
    auto *tree = dynamic_cast<TTree *>(obj);
    if (!tree)
    {
        throw std::runtime_error("Missing entries tree in SampleRootIO directory");
    }

    std::string entry_name;
    std::string artio_path;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot = 0.0;
    double db_tor101_pot = 0.0;
    double normalisation = 1.0;
    double normalised_pot_sum = 0.0;

    tree->SetBranchAddress("entry_name", &entry_name);
    tree->SetBranchAddress("artio_path", &artio_path);
    tree->SetBranchAddress("subrun_pot_sum", &subrun_pot_sum);
    tree->SetBranchAddress("db_tortgt_pot", &db_tortgt_pot);
    tree->SetBranchAddress("db_tor101_pot", &db_tor101_pot);
    tree->SetBranchAddress("normalisation", &normalisation);
    tree->SetBranchAddress("normalised_pot_sum", &normalised_pot_sum);

    const Long64_t n = tree->GetEntries();
    out.entries.reserve(static_cast<size_t>(n));
    for (Long64_t i = 0; i < n; ++i)
    {
        tree->GetEntry(i);
        ProvenanceInput entry;
        entry.entry_name = entry_name;
        entry.artio_path = artio_path;
        entry.subrun_pot_sum = subrun_pot_sum;
        entry.db_tortgt_pot = db_tortgt_pot;
        entry.db_tor101_pot = db_tor101_pot;
        entry.normalisation = normalisation;
        entry.normalised_pot_sum = normalised_pot_sum;
        out.entries.push_back(std::move(entry));
    }

    return out;
}

} // namespace sample

} // namespace nuxsec
