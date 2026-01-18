/**
 *  @file  lib/NuIO/include/NuIO/RunInfoDB.h
 *
 *  @brief SQLite wrapper for run/subrun summary queries
 */

#ifndef NUIO_RUNINFODB_H
#define NUIO_RUNINFODB_H

#include <sqlite3.h>

#include <string>
#include <vector>

#include "NuIO/ArtProvenanceIO.h"

namespace nuio
{

class RunInfoDB
{
  public:
    explicit RunInfoDB(std::string path);
    ~RunInfoDB();

    RunInfoDB(const RunInfoDB &) = delete;
    RunInfoDB &operator=(const RunInfoDB &) = delete;

    RunInfoSums SumRuninfo(const std::vector<RunSubrun> &pairs) const;

  private:
    void Exec(const std::string &sql) const;
    void Prepare(const std::string &sql, sqlite3_stmt **stmt) const;

    std::string db_path_;
    sqlite3 *db_ = nullptr;
};

}

#endif
