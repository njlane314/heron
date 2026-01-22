/* -- C++ -- */
/**
 *  @file  io/include/SubRunInventoryService.hh
 *
 *  @brief Declarations for SubRun inventory helpers.
 */

#ifndef NUXSEC_IO_SUBRUN_INVENTORY_SERVICE_H
#define NUXSEC_IO_SUBRUN_INVENTORY_SERVICE_H

#include <string>
#include <vector>

namespace nuxsec
{

namespace art
{
struct SubrunSummary;
}

/** \brief Build summary information from SubRun trees. */
class SubRunInventoryService
{
  public:
    /// \brief Scan input files for SubRun information.
    static art::SubrunSummary scan_subruns(const std::vector<std::string> &files);
};

}

#endif // NUXSEC_IO_SUBRUN_INVENTORY_SERVICE_H
