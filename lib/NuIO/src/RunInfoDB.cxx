/**
 *  @file  lib/NuIO/src/RunInfoDB.cxx
 *
 *  @brief Implementation of SQLite-backed run info summaries
 */

#include "NuIO/RunInfoDB.h"

#include <stdexcept>
#include <utility>

namespace nuio
{

RunInfoDB::RunInfoDB(std::string path)
    : db_path_(std::move(path))
{
    sqlite3 *db = nullptr;
    const int rc = sqlite3_open_v2(db_path_.c_str(), &db, SQLITE_OPEN_READONLY, nullptr);
    if (rc != SQLITE_OK || !db)
    {
        std::string msg = (db ? sqlite3_errmsg(db) : "sqlite3_open_v2 failed");
        if (db)
            sqlite3_close(db);
        throw std::runtime_error("Failed to open SQLite DB: " + db_path_ + " : " + msg);
    }
    db_ = db;
}

RunInfoDB::~RunInfoDB()
{
    if (db_)
        sqlite3_close(db_);
}

void RunInfoDB::Exec(const std::string &sql) const
{
    char *err = nullptr;
    const int rc = sqlite3_exec(db_, sql.c_str(), nullptr, nullptr, &err);
    if (rc != SQLITE_OK)
    {
        std::string msg = err ? err : sqlite3_errmsg(db_);
        sqlite3_free(err);
        throw std::runtime_error("SQLite exec failed: " + msg);
    }
}

void RunInfoDB::Prepare(const std::string &sql, sqlite3_stmt **stmt) const
{
    const int rc = sqlite3_prepare_v2(db_, sql.c_str(), -1, stmt, nullptr);
    if (rc != SQLITE_OK || !stmt || !(*stmt))
        throw std::runtime_error("SQLite prepare failed: " + std::string(sqlite3_errmsg(db_)));
}

RunInfoSums RunInfoDB::SumRuninfo(const std::vector<RunSubrun> &pairs) const
{
    if (pairs.empty())
        throw std::runtime_error("DB selection is empty (no run/subrun pairs).");

    Exec("CREATE TEMP TABLE IF NOT EXISTS sel(run INTEGER, subrun INTEGER);");
    Exec("DELETE FROM sel;");
    Exec("BEGIN;");

    sqlite3_stmt *ins = nullptr;
    Prepare("INSERT INTO sel(run, subrun) VALUES(?, ?);", &ins);

    for (const auto &p : pairs)
    {
        sqlite3_reset(ins);
        sqlite3_clear_bindings(ins);
        sqlite3_bind_int(ins, 1, p.run);
        sqlite3_bind_int(ins, 2, p.subrun);

        const int rc = sqlite3_step(ins);
        if (rc != SQLITE_DONE)
        {
            const std::string msg = sqlite3_errmsg(db_);
            sqlite3_finalize(ins);
            Exec("ROLLBACK;");
            throw std::runtime_error("SQLite insert failed: " + msg);
        }
    }

    sqlite3_finalize(ins);
    Exec("COMMIT;");

    sqlite3_stmt *q = nullptr;
    Prepare(
        "SELECT "
        "  IFNULL(SUM(r.tortgt), 0.0) AS tortgt_sum, "
        "  IFNULL(SUM(r.tor101), 0.0) AS tor101_sum, "
        "  IFNULL(SUM(r.tor860), 0.0) AS tor860_sum, "
        "  IFNULL(SUM(r.tor875), 0.0) AS tor875_sum, "
        "  IFNULL(SUM(r.EA9CNT), 0)  AS ea9cnt_sum, "
        "  IFNULL(SUM(r.E1DCNT), 0)  AS e1dcnt_sum, "
        "  IFNULL(SUM(r.EXTTrig), 0) AS exttrig_sum, "
        "  IFNULL(SUM(r.Gate1Trig), 0) AS gate1_sum, "
        "  IFNULL(SUM(r.Gate2Trig), 0) AS gate2_sum "
        "FROM runinfo r "
        "JOIN sel USING(run, subrun);",
        &q);

    const int rc = sqlite3_step(q);
    if (rc != SQLITE_ROW)
    {
        const std::string msg = sqlite3_errmsg(db_);
        sqlite3_finalize(q);
        throw std::runtime_error("SQLite sum query failed: " + msg);
    }

    RunInfoSums out{};
    out.n_pairs_loaded = static_cast<long long>(pairs.size());
    out.tortgt_sum = sqlite3_column_double(q, 0);
    out.tor101_sum = sqlite3_column_double(q, 1);
    out.tor860_sum = sqlite3_column_double(q, 2);
    out.tor875_sum = sqlite3_column_double(q, 3);
    out.EA9CNT_sum = sqlite3_column_int64(q, 4);
    out.E1DCNT_sum = sqlite3_column_int64(q, 5);
    out.EXTTrig_sum = sqlite3_column_int64(q, 6);
    out.Gate1Trig_sum = sqlite3_column_int64(q, 7);
    out.Gate2Trig_sum = sqlite3_column_int64(q, 8);

    sqlite3_finalize(q);
    return out;
}

}
