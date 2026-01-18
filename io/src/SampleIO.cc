/* -- C++ -- */
/**
 *  @file  io/src/SampleIO.cc
 *
 *  @brief Implementation for SampleIO aggregation helpers.
 */

#include "SampleIO.hh"

#include <memory>
#include <stdexcept>
#include <utility>

namespace nuxsec
{

Sample SampleIO::aggregate(const std::string &sample_name, const std::vector<std::string> &artio_files)
{
    if (artio_files.empty())
    {
        throw std::runtime_error("Sample aggregation requires at least one ArtIO file.");
    }

    Sample out;
    out.sample_name = sample_name;

    for (const auto &path : artio_files)
    {
        ArtProvenance prov = ArtProvenanceIO::read(path);
        if (out.stages.empty())
        {
            out.kind = prov.kind;
            out.beam = prov.beam;
        }
        else
        {
            if (prov.kind != out.kind)
            {
                throw std::runtime_error("Sample kind mismatch in ArtIO file: " + path);
            }
            if (prov.beam != out.beam)
            {
                throw std::runtime_error("Beam mode mismatch in ArtIO file: " + path);
            }
        }

        SampleStage stage = make_stage(prov, path);
        out.subrun_pot_sum += stage.subrun_pot_sum;
        out.db_tortgt_pot_sum += stage.db_tortgt_pot;
        out.db_tor101_pot_sum += stage.db_tor101_pot;
        out.stages.push_back(std::move(stage));
    }

    out.normalisation = compute_normalisation(out.subrun_pot_sum, out.db_tortgt_pot_sum);
    out.normalised_pot_sum = out.subrun_pot_sum * out.normalisation;

    return out;
}

void SampleIO::write(const Sample &sample, const std::string &out_file)
{
    std::unique_ptr<TFile> f(TFile::Open(out_file.c_str(), "UPDATE"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open merged output file for UPDATE: " + out_file);
    }

    TDirectory *d = f->GetDirectory("SampleIO");
    if (!d)
    {
        d = f->mkdir("SampleIO");
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
        TTree stages("stages", "ArtIO stages included in sample aggregation");

        std::string stage_name;
        std::string artio_path;
        double subrun_pot_sum = 0.0;
        double db_tortgt_pot = 0.0;
        double db_tor101_pot = 0.0;
        double normalisation = 1.0;
        double normalised_pot_sum = 0.0;

        stages.Branch("stage_name", &stage_name);
        stages.Branch("artio_path", &artio_path);
        stages.Branch("subrun_pot_sum", &subrun_pot_sum);
        stages.Branch("db_tortgt_pot", &db_tortgt_pot);
        stages.Branch("db_tor101_pot", &db_tor101_pot);
        stages.Branch("normalisation", &normalisation);
        stages.Branch("normalised_pot_sum", &normalised_pot_sum);

        for (const auto &stage : sample.stages)
        {
            stage_name = stage.stage_name;
            artio_path = stage.artio_path;
            subrun_pot_sum = stage.subrun_pot_sum;
            db_tortgt_pot = stage.db_tortgt_pot;
            db_tor101_pot = stage.db_tor101_pot;
            normalisation = stage.normalisation;
            normalised_pot_sum = stage.normalised_pot_sum;
            stages.Fill();
        }

        stages.Write("stages", TObject::kOverwrite);
    }

    f->Write();
    f->Close();
}

Sample SampleIO::read(const std::string &in_file)
{
    std::unique_ptr<TFile> f(TFile::Open(in_file.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open merged output file for READ: " + in_file);
    }

    TDirectory *d = f->GetDirectory("SampleIO");
    if (!d)
    {
        throw std::runtime_error("Missing SampleIO directory in file: " + in_file);
    }
    d->cd();

    Sample out;
    {
        TObject *obj = d->Get("sample_name");
        auto *named = dynamic_cast<TNamed *>(obj);
        if (!named)
        {
            throw std::runtime_error("Missing sample_name metadata in SampleIO directory");
        }
        out.sample_name = named->GetTitle();
    }

    {
        TObject *obj = d->Get("sample_kind");
        auto *named = dynamic_cast<TNamed *>(obj);
        if (!named)
        {
            throw std::runtime_error("Missing sample_kind metadata in SampleIO directory");
        }
        out.kind = parse_sample_kind(named->GetTitle());
    }

    {
        TObject *obj = d->Get("beam_mode");
        auto *named = dynamic_cast<TNamed *>(obj);
        if (!named)
        {
            throw std::runtime_error("Missing beam_mode metadata in SampleIO directory");
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

    TObject *obj = d->Get("stages");
    auto *tree = dynamic_cast<TTree *>(obj);
    if (!tree)
    {
        throw std::runtime_error("Missing stages tree in SampleIO directory");
    }

    std::string stage_name;
    std::string artio_path;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot = 0.0;
    double db_tor101_pot = 0.0;
    double normalisation = 1.0;
    double normalised_pot_sum = 0.0;

    tree->SetBranchAddress("stage_name", &stage_name);
    tree->SetBranchAddress("artio_path", &artio_path);
    tree->SetBranchAddress("subrun_pot_sum", &subrun_pot_sum);
    tree->SetBranchAddress("db_tortgt_pot", &db_tortgt_pot);
    tree->SetBranchAddress("db_tor101_pot", &db_tor101_pot);
    tree->SetBranchAddress("normalisation", &normalisation);
    tree->SetBranchAddress("normalised_pot_sum", &normalised_pot_sum);

    const Long64_t n = tree->GetEntries();
    out.stages.reserve(static_cast<size_t>(n));
    for (Long64_t i = 0; i < n; ++i)
    {
        tree->GetEntry(i);
        SampleStage stage;
        stage.stage_name = stage_name;
        stage.artio_path = artio_path;
        stage.subrun_pot_sum = subrun_pot_sum;
        stage.db_tortgt_pot = db_tortgt_pot;
        stage.db_tor101_pot = db_tor101_pot;
        stage.normalisation = normalisation;
        stage.normalised_pot_sum = normalised_pot_sum;
        out.stages.push_back(std::move(stage));
    }

    return out;
}

double SampleIO::compute_normalisation(double subrun_pot_sum, double db_tortgt_pot)
{
    if (subrun_pot_sum <= 0.0)
    {
        return 1.0;
    }
    if (db_tortgt_pot <= 0.0)
    {
        return 1.0;
    }
    return db_tortgt_pot / subrun_pot_sum;
}

SampleStage SampleIO::make_stage(const ArtProvenance &prov, const std::string &artio_path)
{
    SampleStage stage;
    stage.stage_name = prov.cfg.stage_name;
    stage.artio_path = artio_path;
    stage.subrun_pot_sum = prov.subrun.pot_sum;
    stage.db_tortgt_pot = prov.db_tortgt_pot;
    stage.db_tor101_pot = prov.db_tor101_pot;
    stage.normalisation = compute_normalisation(stage.subrun_pot_sum, stage.db_tortgt_pot);
    stage.normalised_pot_sum = stage.subrun_pot_sum * stage.normalisation;
    return stage;
}

} // namespace nuxsec
