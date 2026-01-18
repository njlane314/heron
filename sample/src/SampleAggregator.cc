/* -- C++ -- */
/**
 *  @file  sample/src/SampleAggregator.cc
 *
 *  @brief Implementation for SampleAggregator helpers.
 */

#include "SampleAggregator.hh"

#include <memory>
#include <stdexcept>
#include <utility>

namespace nuxsec
{

Sample SampleAggregator::Aggregate(const std::string &sample_name,
                                  const std::vector<std::string> &artio_files)
{
    if (artio_files.empty())
    {
        throw std::runtime_error("Sample aggregation requires at least one art provenance file.");
    }

    Sample out;
    out.sample_name = sample_name;

    for (const auto &path : artio_files)
    {
        ArtFileProvenance prov = ArtFileProvenanceRootIO::read(path);
        if (out.fragments.empty())
        {
            out.kind = prov.kind;
            out.beam = prov.beam;
        }
        else
        {
            if (prov.kind != out.kind)
            {
                throw std::runtime_error("Sample kind mismatch in art provenance file: " + path);
            }
            if (prov.beam != out.beam)
            {
                throw std::runtime_error("Beam mode mismatch in art provenance file: " + path);
            }
        }

        SampleFragment fragment = MakeFragment(prov, path);
        out.subrun_pot_sum += fragment.subrun_pot_sum;
        out.db_tortgt_pot_sum += fragment.db_tortgt_pot;
        out.db_tor101_pot_sum += fragment.db_tor101_pot;
        out.fragments.push_back(std::move(fragment));
    }

    out.normalization = ComputeNormalization(out.subrun_pot_sum, out.db_tortgt_pot_sum);
    out.normalized_pot_sum = out.subrun_pot_sum * out.normalization;

    return out;
}

double SampleAggregator::ComputeNormalization(double subrun_pot_sum, double db_tortgt_pot)
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

SampleFragment SampleAggregator::MakeFragment(const ArtFileProvenance &prov, const std::string &artio_path)
{
    SampleFragment fragment;
    fragment.fragment_name = prov.cfg.stage_name;
    fragment.artio_path = artio_path;
    fragment.subrun_pot_sum = prov.subrun.pot_sum;
    fragment.db_tortgt_pot = prov.db_tortgt_pot;
    fragment.db_tor101_pot = prov.db_tor101_pot;
    fragment.normalization = ComputeNormalization(fragment.subrun_pot_sum, fragment.db_tortgt_pot);
    fragment.normalized_pot_sum = fragment.subrun_pot_sum * fragment.normalization;
    return fragment;
}

} // namespace nuxsec
