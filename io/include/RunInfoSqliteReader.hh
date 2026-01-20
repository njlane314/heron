/* -- C++ -- */
/**
 *  @file  io/include/RunInfoSqliteReader.hh
 *
 *  @brief SQLite reader for run/subrun summary queries.
 */

#ifndef NUXSEC_IO_RUNINFO_SQLITE_READER_H
#define NUXSEC_IO_RUNINFO_SQLITE_READER_H

#include <sqlite3.h>

#include <string>
#include <vector>

#include "ArtFileProvenanceIO.hh"

namespace nuxsec
{

struct RunInfoSums
{
    double tortgt_sum = 0.0;
    double tor101_sum = 0.0;
    double tor860_sum = 0.0;
    double tor875_sum = 0.0;

    long long EA9CNT_sum = 0;
    long long E1DCNT_sum = 0;
    long long EXTTrig_sum = 0;
    long long Gate1Trig_sum = 0;
    long long Gate2Trig_sum = 0;

    long long n_pairs_loaded = 0;
};

class RunInfoSqliteReader
{
  public:
    explicit RunInfoSqliteReader(std::string path);
    ~RunInfoSqliteReader();

    RunInfoSqliteReader(const RunInfoSqliteReader &) = delete;
    RunInfoSqliteReader &operator=(const RunInfoSqliteReader &) = delete;

    RunInfoSums sum_run_info(const std::vector<RunSubrun> &pairs) const;

  private:
    void exec(const std::string &sql) const;
    void prepare(const std::string &sql, sqlite3_stmt **stmt) const;

    std::string db_path_;
    sqlite3 *db_ = nullptr;
};

} // namespace nuxsec

#endif // NUXSEC_IO_RUNINFO_SQLITE_READER_H
