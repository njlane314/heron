/* -- C++ -- */
/**
 *  @file  io/include/SampleTypes.hh
 *
 *  @brief Sample kind and beam mode helpers.
 */

#ifndef Nuxsec_IO_SAMPLE_TYPES_H_INCLUDED
#define Nuxsec_IO_SAMPLE_TYPES_H_INCLUDED

#include <string>

namespace nuxsec
{

enum class SampleKind
{
    kUnknown = 0,
    kData,
    kEXT,
    kOverlay,
    kDirt,
    kStrangeness
};

const char *sample_kind_name(SampleKind k);
SampleKind parse_sample_kind(const std::string &name);

enum class BeamMode
{
    kUnknown = 0,
    kNuMI,
    kBNB
};

const char *beam_mode_name(BeamMode b);
BeamMode parse_beam_mode(const std::string &name);

} // namespace nuxsec

#endif // Nuxsec_IO_SAMPLE_TYPES_H_INCLUDED
