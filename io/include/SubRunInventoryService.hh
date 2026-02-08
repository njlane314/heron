/* -- C++ -- */
/**
 *  @file  io/include/SubRunInventoryService.hh
 *
 *  @brief Declarations for SubRun inventory helpers that track availability,
 *         counts, and metadata summaries.
 */

#ifndef NUXSEC_IO_SUBRUN_INVENTORY_SERVICE_H
#define NUXSEC_IO_SUBRUN_INVENTORY_SERVICE_H

#include <string>
#include <vector>


struct Summary;

class SubRunInventoryService
{
  public:
    static Summary scan_subruns(const std::vector<std::string> &files);
};

#endif // NUXSEC_IO_SUBRUN_INVENTORY_SERVICE_H
