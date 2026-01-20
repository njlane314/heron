/* -- C++ -- */
/**
 *  @file  io/include/RunInfoSqliteReader.hh
 *
 *  @brief SQLite reader for run/subrun summary queries.
 */

#ifndef Nuxsec_IO_RUNINFO_SQLITE_READER_H_INCLUDED
#define Nuxsec_IO_RUNINFO_SQLITE_READER_H_INCLUDED

#include <sqlite3.h>

#include <string>
#include <vector>

#include "ArtFileProvenanceIO.hh"

namespace nuxsec
{

class RunInfoSqliteReader
{
  public:
    explicit RunInfoSqliteReader(std::string path);
    ~RunInfoSqliteReader();

    RunInfoSqliteReader(const RunInfoSqliteReader &) = delete;
    RunInfoSqliteReader &operator=(const RunInfoSqliteReader &) = delete;

    RunInfoSums sumRunInfo(const std::vector<RunSubrun> &pairs) const;

  private:
    void exec(const std::string &sql) const;
    void prepare(const std::string &sql, sqlite3_stmt **stmt) const;

    std::string db_path_;
    sqlite3 *db_ = nullptr;
};

} // namespace nuxsec

#endif // Nuxsec_IO_RUNINFO_SQLITE_READER_H_INCLUDED
