/* -- C++ -- */
/**
 *  @file  io/include/RunInfoDB.hh
 *
 *  @brief SQLite wrapper for run/subrun summary queries.
 */

#ifndef Nuxsec_IO_RUNINFODB_H_INCLUDED
#define Nuxsec_IO_RUNINFODB_H_INCLUDED

#include <sqlite3.h>

#include <string>
#include <vector>

#include "ArtProvenanceIO.hh"

namespace nuxsec
{

class RunInfoDB
{
  public:
    explicit RunInfoDB(std::string path);
    ~RunInfoDB();

    RunInfoDB(const RunInfoDB &) = delete;
    RunInfoDB &operator=(const RunInfoDB &) = delete;

    RunInfoSums sum_runinfo(const std::vector<RunSubrun> &pairs) const;

  private:
    void exec(const std::string &sql) const;
    void prepare(const std::string &sql, sqlite3_stmt **stmt) const;

    std::string db_path_;
    sqlite3 *db_ = nullptr;
};

} // namespace nuxsec

#endif // Nuxsec_IO_RUNINFODB_H_INCLUDED
