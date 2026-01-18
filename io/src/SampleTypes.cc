/* -- C++ -- */
/**
 *  @file  io/src/SampleTypes.cc
 *
 *  @brief Implementation for sample kind and beam mode helpers.
 */

#include "SampleTypes.hh"

#include <algorithm>
#include <cctype>

namespace nuxsec
{

const char *sample_kind_name(SampleKind k)
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

SampleKind parse_sample_kind(const std::string &name)
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

const char *beam_mode_name(BeamMode b)
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

BeamMode parse_beam_mode(const std::string &name)
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

} // namespace nuxsec
