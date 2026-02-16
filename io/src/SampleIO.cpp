/* -- C++ -- */
/**
 *  @file  io/src/SampleIO.cpp
 *
 *  @brief Implementation for SampleIO helpers.
 */

#include "SampleIO.hh"

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <utility>

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include "ArtFileProvenanceIO.hh"



const char *SampleIO::sample_origin_name(SampleOrigin k)
{
    switch (k)
    {
        case SampleOrigin::kData:
            return "data";
        case SampleOrigin::kEXT:
            return "ext";
        case SampleOrigin::kOverlay:
            return "overlay";
        case SampleOrigin::kDirt:
            return "dirt";
        case SampleOrigin::kStrangeness:
            return "strangeness";
        default:
            return "unknown";
    }
}

SampleIO::SampleOrigin SampleIO::parse_sample_origin(const std::string &name)
{
    std::string lowered = name;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                   [](unsigned char c)
                   {
                       return static_cast<char>(std::tolower(c));
                   });

    if (lowered == "data")
    {
        return SampleOrigin::kData;
    }
    if (lowered == "ext")
    {
        return SampleOrigin::kEXT;
    }
    if (lowered == "overlay")
    {
        return SampleOrigin::kOverlay;
    }
    if (lowered == "dirt")
    {
        return SampleOrigin::kDirt;
    }
    if (lowered == "strangeness")
    {
        return SampleOrigin::kStrangeness;
    }
    return SampleOrigin::kUnknown;
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
    TNamed("sample_kind", sample_origin_name(sample.origin)).Write("sample_kind", TObject::kOverwrite);
    TNamed("beam_mode", beam_mode_name(sample.beam)).Write("beam_mode", TObject::kOverwrite);

    TParameter<double>("subrun_pot_sum", sample.subrun_pot_sum).Write("subrun_pot_sum", TObject::kOverwrite);
    TParameter<double>("db_tortgt_pot_sum", sample.db_tortgt_pot_sum).Write("db_tortgt_pot_sum", TObject::kOverwrite);
    TParameter<double>("db_tor101_pot_sum", sample.db_tor101_pot_sum).Write("db_tor101_pot_sum", TObject::kOverwrite);
    TParameter<double>("normalisation", sample.normalisation).Write("normalisation", TObject::kOverwrite);
    TParameter<double>("normalised_pot_sum", sample.normalised_pot_sum).Write("normalised_pot_sum", TObject::kOverwrite);

    {
        TTree entries("entries", "Art file entries included in sample aggregation");

        std::string entry_name;
        std::string art_path;
        double subrun_pot_sum = 0.0;
        double db_tortgt_pot = 0.0;
        double db_tor101_pot = 0.0;
        double normalisation = 1.0;
        double normalised_pot_sum = 0.0;

        entries.Branch("entry_name", &entry_name);
        entries.Branch("art_path", &art_path);
        entries.Branch("subrun_pot_sum", &subrun_pot_sum);
        entries.Branch("db_tortgt_pot", &db_tortgt_pot);
        entries.Branch("db_tor101_pot", &db_tor101_pot);
        entries.Branch("normalisation", &normalisation);
        entries.Branch("normalised_pot_sum", &normalised_pot_sum);

        for (const auto &input : sample.inputs)
        {
            entry_name = input.entry_name;
            art_path = input.art_path;
            subrun_pot_sum = input.subrun_pot_sum;
            db_tortgt_pot = input.db_tortgt_pot;
            db_tor101_pot = input.db_tor101_pot;
            normalisation = input.normalisation;
            normalised_pot_sum = input.normalised_pot_sum;
            entries.Fill();
        }

        entries.Write("entries", TObject::kOverwrite);
    }

    {
        TTree root_files("root_files", "Resolved ROOT input files for sample");

        std::string root_file;
        root_files.Branch("root_file", &root_file);

        for (const auto &path : sample.root_files)
        {
            root_file = path;
            root_files.Fill();
        }

        root_files.Write("root_files", TObject::kOverwrite);
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
        out.origin = parse_sample_origin(named->GetTitle());
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

    std::string *p_entry_name = nullptr;
    std::string *p_art_path = nullptr;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot = 0.0;
    double db_tor101_pot = 0.0;
    double normalisation = 1.0;
    double normalised_pot_sum = 0.0;

    tree->SetBranchAddress("entry_name", &p_entry_name);
    tree->SetBranchAddress("art_path", &p_art_path);
    tree->SetBranchAddress("subrun_pot_sum", &subrun_pot_sum);
    tree->SetBranchAddress("db_tortgt_pot", &db_tortgt_pot);
    tree->SetBranchAddress("db_tor101_pot", &db_tor101_pot);
    tree->SetBranchAddress("normalisation", &normalisation);
    tree->SetBranchAddress("normalised_pot_sum", &normalised_pot_sum);

    const Long64_t n = tree->GetEntries();
    out.inputs.reserve(static_cast<size_t>(n));
    for (Long64_t i = 0; i < n; ++i)
    {
        tree->GetEntry(i);
        if (!p_entry_name || !p_art_path)
        {
            throw std::runtime_error("Missing entry_name or art_path branch data");
        }
        ProvenanceInput input;
        input.entry_name = *p_entry_name;
        input.art_path = *p_art_path;
        input.subrun_pot_sum = subrun_pot_sum;
        input.db_tortgt_pot = db_tortgt_pot;
        input.db_tor101_pot = db_tor101_pot;
        input.normalisation = normalisation;
        input.normalised_pot_sum = normalised_pot_sum;
        out.inputs.push_back(std::move(input));
    }

    TObject *files_obj = d->Get("root_files");
    auto *files_tree = dynamic_cast<TTree *>(files_obj);
    if (files_tree)
    {
        std::string *p_root_file = nullptr;
        files_tree->SetBranchAddress("root_file", &p_root_file);
        const Long64_t n_files = files_tree->GetEntries();
        out.root_files.reserve(static_cast<size_t>(n_files));
        for (Long64_t i = 0; i < n_files; ++i)
        {
            files_tree->GetEntry(i);
            if (!p_root_file)
            {
                throw std::runtime_error("Missing root_file branch data");
            }
            out.root_files.push_back(*p_root_file);
        }
    }

    return out;
}

std::vector<std::string> SampleIO::resolve_root_files(const Sample &sample)
{
    if (!sample.root_files.empty())
    {
        std::vector<std::string> files = sample.root_files;
        std::sort(files.begin(), files.end());
        files.erase(std::unique(files.begin(), files.end()), files.end());
        return files;
    }

    std::vector<std::string> files;
    files.reserve(sample.inputs.size());

    for (const ProvenanceInput &input : sample.inputs)
    {
        const Provenance prov = ArtFileProvenanceIO::read(input.art_path);
        files.insert(files.end(), prov.input_files.begin(), prov.input_files.end());
    }

    std::sort(files.begin(), files.end());
    files.erase(std::unique(files.begin(), files.end()), files.end());

    return files;
}


