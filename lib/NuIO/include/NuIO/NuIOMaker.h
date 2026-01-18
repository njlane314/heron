/**
 *  @file  lib/NuIO/include/NuIO/NuIOMaker.h
 *
 *  @brief Definitions for NuIO maker helpers
 */

#ifndef NUIO_NUIOMAKER_H
#define NUIO_NUIOMAKER_H

#include <string>
#include <vector>

namespace nuio
{

struct NuIOMakerConfig
{
    std::string sample_name;
    std::vector<std::string> artio_files;
    std::string output_path;
};

class NuIOMaker
{
  public:
    static void write(const NuIOMakerConfig &config);
};

} // namespace nuio

#endif
