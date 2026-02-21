/* -- C++ -- */
/**
 *  @file  framework/io/src/NormalisationService.cc
 *
 *  @brief Implementation for Sample normalisation service helpers.
 */

#include "NormalisationService.hh"

#include <stdexcept>
#include <utility>

#include "RunDatabaseService.hh"


SampleIO::Sample NormalisationService::build_sample(const std::string &sample_name,
                                                            const std::vector<std::string> &art_files,
                                                            const std::string &db_path)
{
    if (art_files.empty())
    {
        throw std::runtime_error("Sample aggregation requires at least one Art file provenance root file.");
    }

    SampleIO::Sample out;
    out.sample_name = sample_name;

    RunDatabaseService db(db_path);

    for (const auto &path : art_files)
    {
        Provenance prov = ArtFileProvenanceIO::read(path);
        if (out.inputs.empty())
        {
            out.origin = prov.kind;
            out.beam = prov.beam;
        }
        else
        {
            if (prov.kind != out.origin)
            {
                throw std::runtime_error("Sample kind mismatch in Art file provenance: " + path);
            }
            if (prov.beam != out.beam)
            {
                throw std::runtime_error("Beam mode mismatch in Art file provenance: " + path);
            }
        }

        RunInfoSums runinfo = db.sum_run_info(prov.summary.unique_pairs);
        const double pot_scale = (prov.scale > 0.0) ? prov.scale : 1.0;
        const double db_pot_scale = 1.0e12; // Run DB stores POT in units of 1e12.
        runinfo.tortgt_sum *= pot_scale * db_pot_scale;
        runinfo.tor101_sum *= pot_scale * db_pot_scale;
        runinfo.tor860_sum *= pot_scale * db_pot_scale;
        runinfo.tor875_sum *= pot_scale * db_pot_scale;

        const double db_tortgt_pot = runinfo.tortgt_sum;
        const double db_tor101_pot = runinfo.tor101_sum;

        SampleIO::ProvenanceInput input = make_entry(prov, path, db_tortgt_pot, db_tor101_pot);
        out.subrun_pot_sum += input.subrun_pot_sum;
        out.db_tortgt_pot_sum += input.db_tortgt_pot;
        out.db_tor101_pot_sum += input.db_tor101_pot;
        out.inputs.push_back(std::move(input));
    }

    out.root_files = SampleIO::resolve_root_files(out);
    out.normalisation = compute_normalisation(out.subrun_pot_sum, out.db_tortgt_pot_sum);
    out.normalised_pot_sum = out.subrun_pot_sum * out.normalisation;

    return out;
}

double NormalisationService::compute_normalisation(double subrun_pot_sum, double db_tortgt_pot)
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

SampleIO::ProvenanceInput NormalisationService::make_entry(const Provenance &prov,
                                                                   const std::string &art_path,
                                                                   double db_tortgt_pot,
                                                                   double db_tor101_pot)
{
    SampleIO::ProvenanceInput entry;
    entry.entry_name = prov.input.input_name;
    entry.art_path = art_path;
    entry.subrun_pot_sum = prov.summary.pot_sum;
    entry.db_tortgt_pot = db_tortgt_pot;
    entry.db_tor101_pot = db_tor101_pot;
    entry.normalisation = compute_normalisation(entry.subrun_pot_sum, entry.db_tortgt_pot);
    entry.normalised_pot_sum = entry.subrun_pot_sum * entry.normalisation;
    return entry;
}

