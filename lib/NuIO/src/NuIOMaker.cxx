/**
 *  @file  lib/NuIO/src/NuIOMaker.cxx
 *
 *  @brief Implementation for NuIO maker helpers
 */

#include "NuIO/NuIOMaker.h"

#include <memory>
#include <stdexcept>

#include "TDirectory.h"
#include "TFile.h"
#include "TNamed.h"
#include "TObject.h"
#include "TParameter.h"
#include "TTree.h"

namespace nuio
{

void NuIOMaker::write(const NuIOMakerConfig &config)
{
    if (config.output_path.empty())
    {
        throw std::runtime_error("NuIOMaker output_path is required.");
    }
    if (config.sample_name.empty())
    {
        throw std::runtime_error("NuIOMaker sample_name is required.");
    }
    if (config.artio_files.empty())
    {
        throw std::runtime_error("NuIOMaker requires at least one ArtIO input.");
    }

    std::unique_ptr<TFile> f(TFile::Open(config.output_path.c_str(), "RECREATE"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("Failed to open NuIO output file: " + config.output_path);
    }

    TDirectory *d = f->mkdir("NuIOMaker");
    if (!d)
    {
        throw std::runtime_error("Failed to create NuIOMaker directory in output file.");
    }
    d->cd();

    TNamed("sample_name", config.sample_name.c_str()).Write("sample_name", TObject::kOverwrite);
    TParameter<int>("artio_file_count", static_cast<int>(config.artio_files.size()))
        .Write("artio_file_count", TObject::kOverwrite);

    TTree artio_inputs("artio_inputs", "ArtIO inputs used to assemble NuIO");
    std::string artio_path;
    artio_inputs.Branch("artio_path", &artio_path);
    for (const auto &path : config.artio_files)
    {
        artio_path = path;
        artio_inputs.Fill();
    }
    artio_inputs.Write("artio_inputs", TObject::kOverwrite);

    f->Write();
    f->Close();
}

} // namespace nuio
