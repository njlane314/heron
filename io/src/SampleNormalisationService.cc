/* -- C++ -- */
/**
 *  @file  io/src/SampleNormalisationService.cc
 *
 *  @brief Implementation for Sample normalisation service helpers.
 */

#include "SampleNormalisationService.hh"
#include "RunInfoService.hh"

#include <stdexcept>
#include <utility>

namespace nuxsec
{

SampleIO::Sample SampleNormalisationService::aggregate(const std::string &sample_name,
                                                       const std::vector<std::string> &artio_files,
                                                       const std::string &db_path)
{
    if (artio_files.empty())
    {
        throw std::runtime_error("Sample aggregation requires at least one Art file provenance root file.");
    }

    SampleIO::Sample out;
    out.sample_name = sample_name;

    RunInfoService db(db_path);

    for (const auto &path : artio_files)
    {
        artio::Provenance prov = ArtFileProvenanceIO::read(path);
        if (out.fragments.empty())
        {
            out.kind = prov.kind;
            out.beam = prov.beam;
        }
        else
        {
            if (prov.kind != out.kind)
            {
                throw std::runtime_error("Sample kind mismatch in Art file provenance: " + path);
            }
            if (prov.beam != out.beam)
            {
                throw std::runtime_error("Beam mode mismatch in Art file provenance: " + path);
            }
        }

        RunInfoSums runinfo = db.sum_run_info(prov.subrun.unique_pairs);
        const double pot_scale = (prov.scale > 0.0) ? prov.scale : 1.0;
        runinfo.tortgt_sum *= pot_scale;
        runinfo.tor101_sum *= pot_scale;
        runinfo.tor860_sum *= pot_scale;
        runinfo.tor875_sum *= pot_scale;

        const double db_tortgt_pot = runinfo.tortgt_sum;
        const double db_tor101_pot = runinfo.tor101_sum;

        SampleIO::SampleFragment fragment = make_fragment(prov, path, db_tortgt_pot, db_tor101_pot);
        out.subrun_pot_sum += fragment.subrun_pot_sum;
        out.db_tortgt_pot_sum += fragment.db_tortgt_pot;
        out.db_tor101_pot_sum += fragment.db_tor101_pot;
        out.fragments.push_back(std::move(fragment));
    }

    out.normalisation = compute_normalisation(out.subrun_pot_sum, out.db_tortgt_pot_sum);
    out.normalised_pot_sum = out.subrun_pot_sum * out.normalisation;

    return out;
}

double SampleNormalisationService::compute_normalisation(double subrun_pot_sum, double db_tortgt_pot)
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

SampleIO::SampleFragment SampleNormalisationService::make_fragment(const artio::Provenance &prov,
                                                                   const std::string &artio_path,
                                                                   double db_tortgt_pot,
                                                                   double db_tor101_pot)
{
    SampleIO::SampleFragment fragment;
    fragment.fragment_name = prov.cfg.stage_name;
    fragment.artio_path = artio_path;
    fragment.subrun_pot_sum = prov.subrun.pot_sum;
    fragment.db_tortgt_pot = db_tortgt_pot;
    fragment.db_tor101_pot = db_tor101_pot;
    fragment.normalisation = compute_normalisation(fragment.subrun_pot_sum, fragment.db_tortgt_pot);
    fragment.normalised_pot_sum = fragment.subrun_pot_sum * fragment.normalisation;
    return fragment;
}

} // namespace nuxsec
