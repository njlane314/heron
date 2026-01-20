/* -- C++ -- */
/**
 *  @file  io/include/SubrunTreeScanner.hh
 *
 *  @brief Declarations for Subrun tree scanning helpers.
 */

#ifndef Nuxsec_IO_SUBRUN_TREE_SCANNER_H_INCLUDED
#define Nuxsec_IO_SUBRUN_TREE_SCANNER_H_INCLUDED

#include <string>
#include <vector>

#include "ArtFileProvenanceIO.hh"

namespace nuxsec
{
SubrunTreeSummary scan_subrun_tree(const std::vector<std::string> &files);
}

#endif // Nuxsec_IO_SUBRUN_TREE_SCANNER_H_INCLUDED
