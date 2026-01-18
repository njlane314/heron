/**
 *  @file  lib/NuIO/include/NuIO/ArtIOManifestIO.h
 *
 *  @brief Manifest IO helpers for ArtIO stage metadata
 */

#ifndef NUIO_ARTIO_MANIFEST_IO_H
#define NUIO_ARTIO_MANIFEST_IO_H

#include <string>
#include <vector>

#include "NuIO/ArtProvenanceIO.h"

namespace nuio
{

class ArtIOManifestIO
{
  public:
    static std::vector<std::string> ListStages(const std::string &artio_file);
    static void WriteStage(const std::string &artio_file,
                           const std::string &db_path,
                           double pot_scale,
                           const ArtProvenance &stage);
};

}

#endif
