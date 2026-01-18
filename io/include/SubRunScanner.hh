/* -- C++ -- */
/**
 *  @file  io/include/SubRunScanner.hh
 *
 *  @brief Declarations for SubRun scanning helpers.
 */

#ifndef Nuxsec_IO_SUBRUN_SCANNER_H_INCLUDED
#define Nuxsec_IO_SUBRUN_SCANNER_H_INCLUDED

#include <string>
#include <vector>

#include "ArtProvenanceIO.hh"

namespace nuxsec
{
SubRunInfo scan_subrun_tree(const std::vector<std::string> &files);
}

#endif // Nuxsec_IO_SUBRUN_SCANNER_H_INCLUDED
