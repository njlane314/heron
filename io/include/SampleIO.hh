/* -- C++ -- */
/**
 *  @file  io/include/SampleIO.hh
 *
 *  @brief Sample data structures and ROOT IO helpers.
 */

#ifndef NUXSEC_IO_SAMPLE_IO_H
#define NUXSEC_IO_SAMPLE_IO_H

#include <string>
#include <vector>

namespace nuxsec
{

namespace sample
{

class SampleIO
{
  public:
    enum class SampleKind
    {
        kUnknown = 0,
        kData,
        kEXT,
        kOverlay,
        kDirt,
        kStrangeness
    };

    enum class BeamMode
    {
        kUnknown = 0,
        kNuMI,
        kBNB
    };

    struct ProvenanceInput
    {
        std::string fragment_name;
        std::string artio_path;

        double subrun_pot_sum = 0.0;
        double db_tortgt_pot = 0.0;
        double db_tor101_pot = 0.0;

        double normalisation = 1.0;
        double normalised_pot_sum = 0.0;
    };

    struct Sample
    {
        std::string sample_name;
        SampleKind kind = SampleKind::kUnknown;
        BeamMode beam = BeamMode::kUnknown;

        std::vector<ProvenanceInput> fragments;

        double subrun_pot_sum = 0.0;
        double db_tortgt_pot_sum = 0.0;
        double db_tor101_pot_sum = 0.0;

        double normalisation = 1.0;
        double normalised_pot_sum = 0.0;
    };

    static const char *sample_kind_name(SampleKind k);
    static SampleKind parse_sample_kind(const std::string &name);

    static const char *beam_mode_name(BeamMode b);
    static BeamMode parse_beam_mode(const std::string &name);

    static void write(const Sample &sample, const std::string &out_file);
    static Sample read(const std::string &in_file);
};

} // namespace sample

} // namespace nuxsec

#endif // NUXSEC_IO_SAMPLE_IO_H
