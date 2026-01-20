/* -- C++ -- */
/**
 *  @file  io/include/SubrunTreeScanner.hh
 *
 *  @brief Declarations for Subrun tree scanning helpers.
 */

#ifndef NUXSEC_IO_SUBRUN_TREE_SCANNER_H
#define NUXSEC_IO_SUBRUN_TREE_SCANNER_H

#include <string>
#include <vector>

#include "ArtFileProvenanceIO.hh"

namespace nuxsec
{
SubrunTreeSummary scan_subrun_tree(const std::vector<std::string> &files);
}

#endif // NUXSEC_IO_SUBRUN_TREE_SCANNER_H
