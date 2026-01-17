/**
 *  @file  bin/artIOaggregator/artIOaggregator.cxx
 *
 *  @brief Main implementation for the artIOaggregator executable
 */

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TTree.h>

#include <sqlite3.h>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "NuIO/ArtProvenanceIO.h"

namespace nucond
{

using namespace nuio;

static std::string Trim(std::string s)
{
    auto notspace = [](unsigned char c)
    { return std::isspace(c) == 0; };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
    return s;
}

static void EnsureDirLike(const std::string& path)
{
    (void)path;
}

static std::vector<std::string> ReadFileList(const std::string& filelistPath)
{
    std::ifstream fin(filelistPath);
    if (!fin)
    {
        throw std::runtime_error("Failed to open filelist: " + filelistPath +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    std::vector<std::string> files;
    std::string line;
    while (std::getline(fin, line))
    {
        line = Trim(line);
        if (line.empty())
            continue;
        if (!line.empty() && line[0] == '#')
            continue;
        files.push_back(line);
    }
    if (files.empty())
    {
        throw std::runtime_error("Filelist is empty: " + filelistPath);
    }
    return files;
}

static std::vector<std::string> ReadArgsFile(const std::string& argsPath)
{
    std::ifstream fin(argsPath);
    if (!fin)
    {
        throw std::runtime_error("Failed to open args file: " + argsPath +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    std::vector<std::string> args;
    std::string line;
    while (std::getline(fin, line))
    {
        line = Trim(line);
        if (line.empty())
            continue;
        if (!line.empty() && line[0] == '#')
            continue;
        std::istringstream iss(line);
        std::string tok;
        while (iss >> tok)
        {
            args.push_back(tok);
        }
    }
    return args;
}

static SampleKind InferSampleKind(const std::string& stageName)
{
    const std::string s = ToLower(stageName);
    if (s.find("ext") != std::string::npos)
        return SampleKind::kEXT;
    if (s.find("dirt") != std::string::npos)
        return SampleKind::kMCDirt;
    if (s.find("strange") != std::string::npos)
        return SampleKind::kMCStrangeness;
    if (s.find("beam") != std::string::npos)
        return SampleKind::kMCOverlay;
    if (s.find("data") != std::string::npos)
        return SampleKind::kData;
    return SampleKind::kUnknown;
}

static BeamMode InferBeamMode(const std::string& stageName)
{
    const std::string s = ToLower(stageName);
    if (s.find("numi") != std::string::npos)
        return BeamMode::kNuMI;
    if (s.find("bnb") != std::string::npos)
        return BeamMode::kBNB;
    if (s.find("fhc") != std::string::npos || s.find("rhc") != std::string::npos)
        return BeamMode::kNuMI;
    return BeamMode::kUnknown;
}

struct FilePeek
{
    std::optional<bool> isData;
    std::optional<bool> isNuMI;
};

static FilePeek PeekEventFlags(const std::string& rootFile)
{
    FilePeek out;
    std::unique_ptr<TFile> f(TFile::Open(rootFile.c_str(), "READ"));
    if (!f || f->IsZombie())
        return out;

    TTree* t = dynamic_cast<TTree*>(f->Get("EventSelectionFilter"));
    if (!t)
        return out;

    Bool_t is_data = false;
    Bool_t is_numi = false;

    const bool has_is_data = (t->GetBranch("is_data") != nullptr);
    const bool has_is_numi = (t->GetBranch("is_numi") != nullptr);

    if (!has_is_data && !has_is_numi)
        return out;
    if (t->GetEntries() <= 0)
        return out;

    if (has_is_data)
        t->SetBranchAddress("is_data", &is_data);
    if (has_is_numi)
        t->SetBranchAddress("is_numi", &is_numi);

    t->GetEntry(0);

    if (has_is_data)
        out.isData = static_cast<bool>(is_data);
    if (has_is_numi)
        out.isNuMI = static_cast<bool>(is_numi);
    return out;
}

static uint64_t PackRunSubrun(int run, int subrun)
{
    return (static_cast<uint64_t>(static_cast<uint32_t>(run)) << 32) |
           (static_cast<uint64_t>(static_cast<uint32_t>(subrun)));
}

static SubRunInfo ScanSubRunTree(const std::vector<std::string>& files)
{
    SubRunInfo out;

    const std::vector<std::string> candidates = {"nuselection/SubRun", "SubRun"};
    std::string tree_path;
    std::vector<std::string> missing;
    long long added = 0;

    for (const auto& f : files)
    {
        std::unique_ptr<TFile> file(TFile::Open(f.c_str(), "READ"));
        if (!file || file->IsZombie())
        {
            throw std::runtime_error("Failed to open input ROOT file: " + f);
        }

        for (const auto& name : candidates)
        {
            TTree* tree = dynamic_cast<TTree*>(file->Get(name.c_str()));
            if (tree)
            {
                tree_path = name;
                break;
            }
        }

        if (!tree_path.empty())
            break;
    }

    if (tree_path.empty())
    {
        throw std::runtime_error("No input files contained a SubRun tree.");
    }

    TChain chain(tree_path.c_str());
    for (const auto& f : files)
    {
        std::unique_ptr<TFile> file(TFile::Open(f.c_str(), "READ"));
        if (!file || file->IsZombie())
        {
            throw std::runtime_error("Failed to open input ROOT file: " + f);
        }

        TTree* tree = dynamic_cast<TTree*>(file->Get(tree_path.c_str()));
        if (!tree)
        {
            missing.push_back(f);
            continue;
        }

        if (tree->GetBranch("run") == nullptr || tree->GetBranch("subRun") == nullptr || tree->GetBranch("pot") == nullptr)
        {
            throw std::runtime_error("SubRun tree missing required branches (run, subRun, pot) in file: " + f);
        }

        chain.Add(f.c_str());
        ++added;
    }

    if (!missing.empty())
    {
        std::cerr << "[artIOaggregator] warning: missing SubRun tree in " << missing.size()
                  << " input file(s); skipping.\n";
    }

    if (added == 0)
    {
        throw std::runtime_error("No input files contained a SubRun tree.");
    }

    Int_t run = 0;
    Int_t subRun = 0;
    Double_t pot = 0.0;

    chain.SetBranchAddress("run", &run);
    chain.SetBranchAddress("subRun", &subRun);
    chain.SetBranchAddress("pot", &pot);

    const Long64_t n = chain.GetEntries();
    out.n_entries = static_cast<long long>(n);

    std::unordered_set<uint64_t> seen;
    seen.reserve(static_cast<size_t>(n));

    for (Long64_t i = 0; i < n; ++i)
    {
        chain.GetEntry(i);
        out.pot_sum += static_cast<double>(pot);

        const uint64_t key = PackRunSubrun(run, subRun);
        if (seen.insert(key).second)
        {
            out.unique_pairs.push_back(RunSubrun{static_cast<int>(run), static_cast<int>(subRun)});
        }
    }

    std::sort(out.unique_pairs.begin(), out.unique_pairs.end(),
              [](const RunSubrun& a, const RunSubrun& b)
              {
                  if (a.run != b.run)
                      return a.run < b.run;
                  return a.subrun < b.subrun;
              });

    return out;
}

class RunInfoDB
{
  public:
    explicit RunInfoDB(std::string path)
        : db_path_(std::move(path))
    {
        sqlite3* db = nullptr;
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

    ~RunInfoDB()
    {
        if (db_)
            sqlite3_close(db_);
    }

    RunInfoDB(const RunInfoDB&) = delete;
    RunInfoDB& operator=(const RunInfoDB&) = delete;

    RunInfoSums SumRuninfoForSelection(const std::vector<RunSubrun>& pairs) const
    {
        if (pairs.empty())
            throw std::runtime_error("DB selection is empty (no run/subrun pairs).");

        Exec("CREATE TEMP TABLE IF NOT EXISTS sel(run INTEGER, subrun INTEGER);");
        Exec("DELETE FROM sel;");

        Exec("BEGIN;");
        sqlite3_stmt* ins = nullptr;
        Prepare("INSERT INTO sel(run, subrun) VALUES(?, ?);", &ins);

        for (const auto& p : pairs)
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

        RunInfoSums out{};
        out.n_pairs_loaded = static_cast<long long>(pairs.size());

        sqlite3_stmt* q = nullptr;
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

  private:
    void Exec(const std::string& sql) const
    {
        char* err = nullptr;
        const int rc = sqlite3_exec(db_, sql.c_str(), nullptr, nullptr, &err);
        if (rc != SQLITE_OK)
        {
            std::string msg = err ? err : sqlite3_errmsg(db_);
            sqlite3_free(err);
            throw std::runtime_error("SQLite exec failed: " + msg + " ; SQL=" + sql);
        }
    }

    void Prepare(const std::string& sql, sqlite3_stmt** stmt) const
    {
        const int rc = sqlite3_prepare_v2(db_, sql.c_str(), -1, stmt, nullptr);
        if (rc != SQLITE_OK || !stmt || !(*stmt))
        {
            std::string msg = sqlite3_errmsg(db_);
            throw std::runtime_error("SQLite prepare failed: " + msg + " ; SQL=" + sql);
        }
    }

    std::string db_path_;
    sqlite3* db_ = nullptr;
};

struct CLI
{
    std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
    std::string outdir = "./condensed";
    std::string manifest_path = "./condensed/manifest.root";
    double pot_scale = 1e12;
    std::string ext_denom = "EXTTrig";

    std::vector<StageCfg> stages;
};

static void PrintUsage(const char* argv0)
{
    std::cerr << "Usage:\n"
                 "  "
              << argv0 << " [options] --stage NAME:FILELIST [--stage NAME:FILELIST ...]\n\n"
                          "Options:\n"
                          "  --args FILE         Read additional CLI arguments from FILE\n"
                          "  --db PATH           Path to run.db (default: /exp/uboone/data/uboonebeam/beamdb/run.db)\n"
                          "  --outdir DIR        Output directory (default: ./condensed)\n"
                          "  --manifest PATH     Manifest ROOT file (default: ./condensed/manifest.root)\n"
                          "  --pot-scale X       Multiply tortgt/tor101 sums by X to get POT (default: 1e12)\n"
                          "  --ext-denom COL     EXT denominator column (EXTTrig|Gate1Trig|Gate2Trig...) (default: EXTTrig)\n"
                          "  --help              Print this message\n";
}

static std::vector<std::string> ExpandArgs(int argc, char** argv)
{
    std::vector<std::string> expanded;
    expanded.reserve(static_cast<size_t>(argc));
    for (int i = 0; i < argc; ++i)
    {
        if (i > 0 && std::string(argv[i]) == "--args")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("Missing value for --args");
            }
            const std::vector<std::string> fileArgs = ReadArgsFile(argv[++i]);
            expanded.insert(expanded.end(), fileArgs.begin(), fileArgs.end());
            continue;
        }
        expanded.emplace_back(argv[i]);
    }
    return expanded;
}

static CLI ParseArgs(const std::vector<std::string>& args)
{
    CLI cli;
    for (size_t i = 1; i < args.size(); ++i)
    {
        const std::string& a = args[i];
        auto need = [&](const std::string& opt) -> std::string
        {
            if (i + 1 >= args.size())
                throw std::runtime_error("Missing value for " + opt);
            return args[++i];
        };

        if (a == "--help" || a == "-h")
        {
            PrintUsage(args[0].c_str());
            std::exit(0);
        }
        else if (a == "--args")
        {
            throw std::runtime_error("Unexpected --args: use --args only on the command line.");
        }
        else if (a == "--db")
        {
            cli.db_path = need("--db");
        }
        else if (a == "--outdir")
        {
            cli.outdir = need("--outdir");
        }
        else if (a == "--manifest")
        {
            cli.manifest_path = need("--manifest");
        }
        else if (a == "--pot-scale")
        {
            cli.pot_scale = std::stod(need("--pot-scale"));
        }
        else if (a == "--ext-denom")
        {
            cli.ext_denom = need("--ext-denom");
        }
        else if (a == "--stage")
        {
            const std::string spec = need("--stage");
            const auto pos = spec.find(':');
            if (pos == std::string::npos)
            {
                throw std::runtime_error("Bad --stage spec (expected NAME:FILELIST): " + spec);
            }
            StageCfg sc;
            sc.stage_name = spec.substr(0, pos);
            sc.filelist_path = spec.substr(pos + 1);
            sc.stage_name = Trim(sc.stage_name);
            sc.filelist_path = Trim(sc.filelist_path);
            if (sc.stage_name.empty() || sc.filelist_path.empty())
            {
                throw std::runtime_error("Bad --stage spec (empty name or filelist): " + spec);
            }
            cli.stages.push_back(sc);
        }
        else
        {
            throw std::runtime_error("Unknown argument: " + a);
        }
    }

    if (cli.stages.empty())
    {
        throw std::runtime_error("No --stage specified.");
    }
    return cli;
}

static long long GetDenomFromDB(const RunInfoSums& s, const std::string& denomCol)
{
    const std::string c = ToLower(denomCol);
    if (c == "exttrig")
        return s.EXTTrig_sum;
    if (c == "gate1trig")
        return s.Gate1Trig_sum;
    if (c == "gate2trig")
        return s.Gate2Trig_sum;
    return s.EXTTrig_sum;
}

static ArtProvenance ProcessStage(const StageCfg& cfg,
                                  const RunInfoDB& db,
                                  double pot_scale,
                                  const std::string& ext_denom_col,
                                  const std::string& outdir)
{
    ArtProvenance r;
    r.cfg = cfg;

    r.input_files = ReadFileList(cfg.filelist_path);

    r.kind = InferSampleKind(cfg.stage_name);
    r.beam = InferBeamMode(cfg.stage_name);

    const FilePeek peek = PeekEventFlags(r.input_files.front());
    if (peek.isNuMI.has_value())
    {
        r.beam = (*peek.isNuMI ? BeamMode::kNuMI : BeamMode::kBNB);
    }
    if (peek.isData.has_value())
    {
        if (*peek.isData)
        {
            if (r.kind == SampleKind::kUnknown)
                r.kind = SampleKind::kData;
        }
        else
        {
            if (r.kind == SampleKind::kUnknown)
                r.kind = SampleKind::kMCOverlay;
        }
    }

    r.subrun = ScanSubRunTree(r.input_files);

    r.runinfo = db.SumRuninfoForSelection(r.subrun.unique_pairs);

    r.db_tortgt_pot = r.runinfo.tortgt_sum * pot_scale;
    r.db_tor101_pot = r.runinfo.tor101_sum * pot_scale;

    r.scale = 1.0;

    if (r.kind == SampleKind::kMCOverlay || r.kind == SampleKind::kMCDirt || r.kind == SampleKind::kMCStrangeness)
    {
        const double pot_mc = r.subrun.pot_sum;
        const double pot_data = r.db_tortgt_pot;

        if (pot_mc > 0.0 && pot_data > 0.0)
        {
            r.scale = pot_data / pot_mc;
        }
        else
        {
            r.scale = 0.0;
        }
    }
    else if (r.kind == SampleKind::kEXT)
    {
        const long long num = r.runinfo.EA9CNT_sum;
        const long long den = GetDenomFromDB(r.runinfo, ext_denom_col);
        if (num > 0 && den > 0)
            r.scale = static_cast<double>(num) / static_cast<double>(den);
        else
            r.scale = 0.0;
    }
    else
    {
        r.scale = 1.0;
    }

    const std::string outFile = outdir + "/" + cfg.stage_name + ".condensed.root";

    ArtProvenanceIO::Write(r, outFile);

    return r;
}

} // namespace nucond

int main(int argc, char** argv)
{
    using namespace nucond;

    try
    {
        const std::vector<std::string> expandedArgs = ExpandArgs(argc, argv);
        const CLI cli = ParseArgs(expandedArgs);

        EnsureDirLike(cli.outdir);

        RunInfoDB db(cli.db_path);

        std::vector<ArtProvenance> results;
        results.reserve(cli.stages.size());

        for (const auto& st : cli.stages)
        {
            std::cerr << "[artIOaggregator] stage=" << st.stage_name
                      << " filelist=" << st.filelist_path << "\n";
            ArtProvenance r = ProcessStage(st, db, cli.pot_scale, cli.ext_denom, cli.outdir);

            std::cerr << "  kind=" << SampleKindName(r.kind)
                      << " beam=" << BeamModeName(r.beam)
                      << " files=" << r.input_files.size()
                      << " unique_pairs=" << r.subrun.unique_pairs.size()
                      << " subrun_pot_sum=" << r.subrun.pot_sum
                      << " db_tortgt_pot=" << r.db_tortgt_pot
                      << " EA9CNT=" << r.runinfo.EA9CNT_sum
                      << " EXTTrig=" << r.runinfo.EXTTrig_sum
                      << " scale=" << r.scale
                      << "\n";

            results.push_back(std::move(r));
        }

        {
            std::unique_ptr<TFile> mf(TFile::Open(cli.manifest_path.c_str(), "RECREATE"));
            if (!mf || mf->IsZombie())
            {
                throw std::runtime_error("Failed to create manifest file: " + cli.manifest_path);
            }

            TTree t("Stages", "Condenser stage summary");
            std::string stage_name;
            std::string kind;
            std::string beam;
            std::string condensed_file;
            double subrun_pot_sum = 0.0;
            double db_tortgt_pot = 0.0;
            long long ea9cnt_sum = 0;
            long long exttrig_sum = 0;
            double scale = 1.0;
            long long n_files = 0;
            long long n_pairs = 0;

            t.Branch("stage_name", &stage_name);
            t.Branch("kind", &kind);
            t.Branch("beam", &beam);
            t.Branch("condensed_file", &condensed_file);
            t.Branch("subrun_pot_sum", &subrun_pot_sum);
            t.Branch("db_tortgt_pot", &db_tortgt_pot);
            t.Branch("ea9cnt_sum", &ea9cnt_sum);
            t.Branch("exttrig_sum", &exttrig_sum);
            t.Branch("scale", &scale);
            t.Branch("n_files", &n_files);
            t.Branch("n_pairs", &n_pairs);

            for (const auto& r : results)
            {
                stage_name = r.cfg.stage_name;
                kind = SampleKindName(r.kind);
                beam = BeamModeName(r.beam);
                condensed_file = cli.outdir + "/" + r.cfg.stage_name + ".condensed.root";
                subrun_pot_sum = r.subrun.pot_sum;
                db_tortgt_pot = r.db_tortgt_pot;
                ea9cnt_sum = r.runinfo.EA9CNT_sum;
                exttrig_sum = r.runinfo.EXTTrig_sum;
                scale = r.scale;
                n_files = static_cast<long long>(r.input_files.size());
                n_pairs = static_cast<long long>(r.subrun.unique_pairs.size());
                t.Fill();
            }

            mf->cd();
            t.Write();

            TDirectory* d = mf->mkdir("NuCondenser");
            d->cd();
            TNamed("db_path", cli.db_path.c_str()).Write("db_path", TObject::kOverwrite);
            TParameter<double>("pot_scale", cli.pot_scale).Write("pot_scale", TObject::kOverwrite);
            TNamed("ext_denom_column", cli.ext_denom.c_str()).Write("ext_denom_column", TObject::kOverwrite);

            mf->Write();
            mf->Close();
        }

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
