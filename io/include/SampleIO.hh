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
    enum class SampleOrigin
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
        std::string entry_name;
        std::string art_path;

        double subrun_pot_sum = 0.0;
        double db_tortgt_pot = 0.0;
        double db_tor101_pot = 0.0;

        double normalisation = 1.0;
        double normalised_pot_sum = 0.0;
    };

    struct Sample
    {
        std::string sample_name;
        SampleOrigin origin = SampleOrigin::kUnknown;
        BeamMode beam = BeamMode::kUnknown;

        std::vector<ProvenanceInput> inputs;

        double subrun_pot_sum = 0.0;
        double db_tortgt_pot_sum = 0.0;
        double db_tor101_pot_sum = 0.0;

        double normalisation = 1.0;
        double normalised_pot_sum = 0.0;
    };

    static const char *sample_origin_name(SampleOrigin k);
    static SampleOrigin parse_sample_origin(const std::string &name);

    static const char *beam_mode_name(BeamMode b);
    static BeamMode parse_beam_mode(const std::string &name);

    static void write(const Sample &sample, const std::string &out_file);
    static Sample read(const std::string &in_file);
};

} // namespace sample

} // namespace nuxsec

#endif // NUXSEC_IO_SAMPLE_IO_H
