/**
 *  @file  lib/NuIO/include/NuIO/SubRunScanner.h
 *
 *  @brief Declarations for SubRun scanning helpers
 */

#ifndef NUIO_SUBRUN_SCANNER_H
#define NUIO_SUBRUN_SCANNER_H

#include <string>
#include <vector>

#include "NuIO/ArtProvenanceIO.h"

namespace nuio
{
SubRunInfo scan_subrun_tree(const std::vector<std::string> &files);
}

#endif
