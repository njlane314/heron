/* -- C++ -- */
/**
 *  @file  apps/src/nuxsec.cc
 *
 *  @brief Unified CLI for Nuxsec utilities.
 */

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>

#include "AppUtils.hh"
#include "ArtCLI.hh"
#include "EventCLI.hh"
#include "SampleCLI.hh"

namespace nuxsec
{

namespace app
{

namespace art
{

int run(const Args &art_args, const std::string &log_prefix)
{
    const double pot_scale = 1e12;

    std::filesystem::path out_path(art_args.art_path);
    if (!out_path.parent_path().empty())
    {
        std::filesystem::create_directories(out_path.parent_path());
    }

    const auto files = nuxsec::app::read_paths(art_args.input.filelist_path);

    nuxsec::art::Provenance rec;
    rec.input = art_args.input;
    rec.input_files = files;
    rec.kind = art_args.sample_origin;
    rec.beam = art_args.beam_mode;

    if (rec.kind == nuxsec::sample::SampleIO::SampleOrigin::kUnknown &&
        is_selection_data_file(files.front()))
    {
        rec.kind = nuxsec::sample::SampleIO::SampleOrigin::kData;
    }

    const auto start_time = std::chrono::steady_clock::now();
    log_scan_start(log_prefix);
    rec.summary = nuxsec::SubRunInventoryService::scan_subruns(files);
    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    log_scan_finish(log_prefix, rec.summary.n_entries, elapsed_seconds);

    rec.summary.pot_sum *= pot_scale;
    rec.scale = pot_scale;

    std::cerr << "[" << log_prefix << "] add input=" << rec.input.input_name
              << " files=" << rec.input_files.size()
              << " pairs=" << rec.summary.unique_pairs.size()
              << " pot_sum=" << rec.summary.pot_sum
              << "\n";

    nuxsec::ArtFileProvenanceIO::write(rec, art_args.art_path);

    return 0;
}

} // namespace art

namespace sample
{

int run(const Args &sample_args, const std::string &log_prefix)
{
    const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
    const auto files = nuxsec::app::read_paths(sample_args.filelist_path);

    std::filesystem::path output_path(sample_args.output_path);
    if (!output_path.parent_path().empty())
    {
        std::filesystem::create_directories(output_path.parent_path());
    }
    std::filesystem::path sample_list_path(sample_args.sample_list_path);
    if (!sample_list_path.parent_path().empty())
    {
        std::filesystem::create_directories(sample_list_path.parent_path());
    }

    const auto start_time = std::chrono::steady_clock::now();
    log_sample_start(log_prefix, files.size());
    nuxsec::sample::SampleIO::Sample sample =
        nuxsec::NormalisationService::build_sample(sample_args.sample_name, files, db_path);
    sample.root_files = nuxsec::sample::SampleIO::resolve_root_files(sample);
    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    log_sample_finish(log_prefix, sample.inputs.size(), elapsed_seconds);
    nuxsec::sample::SampleIO::write(sample, sample_args.output_path);
    update_sample_list(sample_args.sample_list_path, sample, sample_args.output_path);

    std::cerr << "[" << log_prefix << "] sample=" << sample.sample_name
              << " inputs=" << sample.inputs.size()
              << " pot_sum=" << sample.subrun_pot_sum
              << " db_tortgt_pot_sum=" << sample.db_tortgt_pot_sum
              << " normalisation=" << sample.normalisation
              << " normalised_pot_sum=" << sample.normalised_pot_sum
              << " output=" << sample_args.output_path
              << " sample_list=" << sample_args.sample_list_path
              << "\n";

    return 0;
}

} // namespace sample

namespace event
{

int run(const Args &event_args, const std::string &log_prefix)
{
    if (event_args.nthreads > 0)
        ROOT::EnableImplicitMT(event_args.nthreads);
    else
        ROOT::EnableImplicitMT();

    const auto &analysis = nuxsec::AnalysisConfigService::instance();
    const auto entries = nuxsec::app::read_samples(event_args.list_path);
    const auto start_time = std::chrono::steady_clock::now();
    log_event_start(log_prefix, entries.size());

    std::vector<Input> inputs;
    inputs.reserve(entries.size());
    std::vector<nuxsec::event::SampleInfo> sample_refs;
    sample_refs.reserve(entries.size());

    for (const auto &entry : entries)
    {
        nuxsec::sample::SampleIO::Sample sample = nuxsec::sample::SampleIO::read(entry.output_path);

        nuxsec::event::SampleInfo ref;
        ref.sample_name = sample.sample_name;
        ref.sample_rootio_path = entry.output_path;
        ref.sample_origin = static_cast<int>(sample.origin);
        ref.beam_mode = static_cast<int>(sample.beam);
        ref.subrun_pot_sum = sample.subrun_pot_sum;
        ref.db_tortgt_pot_sum = sample.db_tortgt_pot_sum;
        ref.db_tor101_pot_sum = sample.db_tor101_pot_sum;
        sample_refs.push_back(std::move(ref));

        Input input;
        input.entry = entry;
        input.sample = std::move(sample);
        inputs.push_back(std::move(input));
    }

    const std::string x_column = "evt_x";
    const std::string y_column;
    const std::vector<std::string> double_columns = {
        "w_nominal",
        "w_template",
        "reco_nu_energy"
    };
    const std::vector<std::string> int_columns = {
        "run_i",
        "subrun_i",
        "event_i",
        "evt_cutmask_i",
        "evt_category",
        "sel_template_i",
        "sel_reco_fv_i",
        "sel_signal_i",
        "sel_bkg_i",
        "in_reco_fiducial_i",
        "analysis_channels",
        "scattering_mode",
        "is_strange_i",
        "count_strange",
        "is_signal_i",
        "recognised_signal_i"
    };

    std::ostringstream schema;
    schema << "type\tname\n";
    schema << "x\t" << x_column << "\n";
    if (!y_column.empty())
        schema << "y\t" << y_column << "\n";
    for (const auto &c : double_columns)
        schema << "double\t" << c << "\n";
    for (const auto &c : int_columns)
        schema << "int\t" << c << "\n";

    const std::string provenance_tree = "nuxsec_art_provenance/run_subrun";
    const std::string event_tree = analysis.tree_name();

    nuxsec::event::Header header;
    header.analysis_name = analysis.name();
    header.provenance_tree = provenance_tree;
    header.event_tree = event_tree;
    header.sample_list_source = event_args.list_path;

    const std::filesystem::path output_path(event_args.output_root);
    if (!output_path.parent_path().empty())
        std::filesystem::create_directories(output_path.parent_path());

    nuxsec::event::EventIO::init(event_args.output_root, header, sample_refs, schema.str(), "compiled");
    nuxsec::event::EventIO event_io(event_args.output_root, nuxsec::event::EventIO::OpenMode::kUpdate);

    for (const auto &input : inputs)
    {
        const nuxsec::sample::SampleIO::Sample &sample = input.sample;
        std::cerr << "[" << log_prefix << "]"
                  << " stage=ensure_tree sample=" << sample.sample_name
                  << " tree=" << event_tree
                  << "\n";
        ensure_tree_present(sample, event_tree);
        std::cerr << "[" << log_prefix << "]"
                  << " stage=load_rdf sample=" << sample.sample_name
                  << "\n";
        ROOT::RDataFrame rdf = nuxsec::RDataFrameFactory::load_sample(sample, event_tree);
        std::cerr << "[" << log_prefix << "]"
                  << " stage=make_processor sample=" << sample.sample_name
                  << "\n";
        const nuxsec::ProcessorEntry proc_entry = analysis.make_processor(sample);

        const auto &processor = nuxsec::ColumnDerivationService::instance();
        std::cerr << "[" << log_prefix << "]"
                  << " stage=define_columns sample=" << sample.sample_name
                  << "\n";
        ROOT::RDF::RNode node = processor.define(rdf, proc_entry);

        {
            const auto cnames = node.GetColumnNames();
            auto has = [&](const std::string &name) {
                return std::find(cnames.begin(), cnames.end(), name) != cnames.end();
            };
            if (!has("run"))
                node = node.Define("run", [] { return static_cast<int>(0); });
            if (!has("subrun"))
                node = node.Define("subrun", [] { return static_cast<int>(0); });
            if (!has("event"))
                node = node.Define("event", [] { return static_cast<int>(0); });

            if (!has("reco_nu_energy"))
                node = node.Define("reco_nu_energy", [] { return 0.0f; });
        }

        node = node.Define(
            "evt_cutmask",
            [](bool t, bool fv, bool sig) {
                uint64_t m = 0;
                if (t)   m |= (1ull << 0);
                if (fv)  m |= (1ull << 1);
                if (sig) m |= (1ull << 2);
                return m;
            },
            {"sel_template", "in_reco_fiducial", "recognised_signal"});

        node = node.Define(
            "evt_category",
            [](bool s, bool b) { return s ? 1 : (b ? 0 : -1); },
            {"sel_signal", "sel_bkg"});

        node = node.Define("evt_x", [](float e) { return static_cast<double>(e); }, {"reco_nu_energy"});
        node = node.Define("run_i", "static_cast<int>(run)");
        node = node.Define("subrun_i", "static_cast<int>(subrun)");
        node = node.Define("event_i", "static_cast<int>(event)");
        node = node.Define("evt_cutmask_i", [](uint64_t v) { return static_cast<int>(v); }, {"evt_cutmask"});
        node = node.Define("sel_template_i", [](bool v) { return static_cast<int>(v); }, {"sel_template"});
        node = node.Define("sel_reco_fv_i", [](bool v) { return static_cast<int>(v); }, {"sel_reco_fv"});
        node = node.Define("sel_signal_i", [](bool v) { return static_cast<int>(v); }, {"sel_signal"});
        node = node.Define("sel_bkg_i", [](bool v) { return static_cast<int>(v); }, {"sel_bkg"});
        node = node.Define("in_reco_fiducial_i", [](bool v) { return static_cast<int>(v); }, {"in_reco_fiducial"});
        node = node.Define("is_strange_i", [](bool v) { return static_cast<int>(v); }, {"is_strange"});
        node = node.Define("is_signal_i", [](bool v) { return static_cast<int>(v); }, {"is_signal"});
        node = node.Define("recognised_signal_i", [](bool v) { return static_cast<int>(v); }, {"recognised_signal"});

        using SampleOrigin = nuxsec::sample::SampleIO::SampleOrigin;
        const auto origin = sample.origin;
        if (origin == SampleOrigin::kOverlay)
        {
            std::cerr << "[" << log_prefix << "]"
                      << " stage=filter_overlay sample=" << sample.sample_name
                      << "\n";
            node = node.Filter([](int strange) { return strange == 0; }, {"count_strange"});
        }
        else if (origin == SampleOrigin::kStrangeness)
        {
            std::cerr << "[" << log_prefix << "]"
                      << " stage=filter_strangeness sample=" << sample.sample_name
                      << "\n";
            node = node.Filter([](int strange) { return strange > 0; }, {"count_strange"});
        }

        std::cerr << "[" << log_prefix << "]"
                  << " stage=snapshot sample=" << sample.sample_name
                  << "\n";
        const ULong64_t n_written =
            event_io.snapshot_event_list(node,
                                         sample.sample_name,
                                         x_column,
                                         y_column,
                                         double_columns,
                                         int_columns,
                                         "sel_template");

        std::cerr << "[" << log_prefix << "] analysis=" << analysis.name()
                  << " sample=" << sample.sample_name
                  << " kind=" << nuxsec::sample::SampleIO::sample_origin_name(sample.origin)
                  << " beam=" << nuxsec::sample::SampleIO::beam_mode_name(sample.beam)
                  << " events_written=" << n_written
                  << " output=" << event_args.output_root
                  << "\n";
    }
    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    log_event_finish(log_prefix, entries.size(), elapsed_seconds);

    return 0;
}

} // namespace event

} // namespace app

} // namespace nuxsec

namespace
{

const char *kUsageArt = "Usage: nuxsec art INPUT_NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]";
const char *kUsageEvent = "Usage: nuxsec event SAMPLE_LIST.tsv OUTPUT.root [NTHREADS]";
const char *kUsageSample = "Usage: nuxsec sample NAME:FILELIST";
const char *kUsageMacro =
    "Usage: nuxsec macro MACRO.C [CALL]\n"
    "       nuxsec macro list\n"
    "\nEnvironment:\n"
    "  NUXSEC_PLOT_DIR     Output directory (default: <repo>/build/plot)\n"
    "  NUXSEC_PLOT_FORMAT  Output extension (default: pdf)\n";

bool is_help_arg(const std::string &arg)
{
    return arg == "-h" || arg == "--help";
}

void print_main_help(std::ostream &out)
{
    out << "Usage: nuxsec <command> [args]\n\n"
        << "Commands:\n"
        << "  art         Aggregate art provenance for an input\n"
        << "  sample      Aggregate Sample ROOT files from art provenance\n"
        << "  event       Build event-level output from aggregated samples\n"
        << "  macro       Run plot macros\n"
        << "\nRun 'nuxsec <command> --help' for command-specific usage.\n";
}

std::filesystem::path find_repo_root()
{
    std::vector<std::filesystem::path> candidates;

    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec)
    {
        candidates.push_back(exe.parent_path());
    }
    candidates.push_back(std::filesystem::current_path());

    for (auto base : candidates)
    {
        for (int i = 0; i < 6; ++i)
        {
            if (std::filesystem::exists(base / "plot/macro/.plot_driver.retired"))
            {
                return base;
            }
            if (!base.has_parent_path())
            {
                break;
            }
            base = base.parent_path();
        }
    }

    return std::filesystem::current_path();
}

void ensure_plot_env(const std::filesystem::path &repo_root)
{
    if (!gSystem->Getenv("NUXSEC_REPO_ROOT"))
    {
        gSystem->Setenv("NUXSEC_REPO_ROOT", repo_root.string().c_str());
    }
    if (!gSystem->Getenv("NUXSEC_PLOT_DIR"))
    {
        const auto out = (repo_root / "build" / "plot").string();
        gSystem->Setenv("NUXSEC_PLOT_DIR", out.c_str());
    }
}

int run_art_command(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageArt << "\n";
        return 0;
    }

    const nuxsec::app::art::Args art_args =
        nuxsec::app::art::parse_args(args, kUsageArt);
    return nuxsec::app::art::run(art_args, "nuxsec:art");
}

int run_sample_command(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageSample << "\n";
        return 0;
    }

    const nuxsec::app::sample::Args sample_args =
        nuxsec::app::sample::parse_args(args, kUsageSample);
    return nuxsec::app::sample::run(sample_args, "nuxsec sample");
}

int run_event_command(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageEvent << "\n";
        return 0;
    }

    const nuxsec::app::event::Args event_args =
        nuxsec::app::event::parse_args(args, kUsageEvent);
    return nuxsec::app::event::run(event_args, "nuxsec event");
}


std::filesystem::path resolve_macro_path(const std::filesystem::path &repo_root,
                                         const std::string &macro_path)
{
    std::filesystem::path candidate(macro_path);
    if (candidate.is_relative())
    {
        const auto repo_candidate = repo_root / candidate;
        if (std::filesystem::exists(repo_candidate))
        {
            return repo_candidate;
        }
        const auto macro_candidate = repo_root / "plot" / "macro" / candidate;
        if (std::filesystem::exists(macro_candidate))
        {
            return macro_candidate;
        }
    }
    return candidate;
}

void add_plot_include_paths(const std::filesystem::path &repo_root)
{
    const auto include_path = repo_root / "plot/include";
    gSystem->AddIncludePath(("-I" + include_path.string()).c_str());
    const auto ana_include_path = repo_root / "ana/include";
    gSystem->AddIncludePath(("-I" + ana_include_path.string()).c_str());
}

int run_root_macro_call(const std::filesystem::path &repo_root,
                        const std::filesystem::path &macro_path,
                        const std::string &call_cmd)
{
    ensure_plot_env(repo_root);
    add_plot_include_paths(repo_root);

    if (!std::filesystem::exists(macro_path))
    {
        throw std::runtime_error("Macro not found at " + macro_path.string());
    }

    const std::string load_cmd = ".L " + macro_path.string();
    gROOT->ProcessLine(load_cmd.c_str());
    const long result = call_cmd.empty() ? 0 : gROOT->ProcessLine(call_cmd.c_str());
    return static_cast<int>(result);
}

int run_root_macro_exec(const std::filesystem::path &repo_root,
                        const std::filesystem::path &macro_path)
{
    ensure_plot_env(repo_root);
    add_plot_include_paths(repo_root);

    if (!std::filesystem::exists(macro_path))
    {
        throw std::runtime_error("Macro not found at " + macro_path.string());
    }

    const std::string exec_cmd = ".x " + macro_path.string();
    const long result = gROOT->ProcessLine(exec_cmd.c_str());
    return static_cast<int>(result);
}

void print_macro_list(std::ostream &out, const std::filesystem::path &repo_root)
{
    const auto macro_dir = repo_root / "plot" / "macro";
    out << "Plot macros in " << macro_dir.string() << ":\n";
    if (!std::filesystem::exists(macro_dir))
    {
        out << "  (none; directory not found)\n";
        return;
    }

    std::vector<std::string> macros;
    for (const auto &entry : std::filesystem::directory_iterator(macro_dir))
    {
        if (!entry.is_regular_file())
        {
            continue;
        }
        const auto &path = entry.path();
        if (path.extension() == ".C")
        {
            macros.push_back(path.filename().string());
        }
    }

    std::sort(macros.begin(), macros.end());
    for (const auto &macro : macros)
    {
        out << "  " << macro << "\n";
    }
}

int run_macro_command(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageMacro << "\n";
        print_macro_list(std::cout, find_repo_root());
        return 0;
    }

    const auto repo_root = find_repo_root();
    ensure_plot_env(repo_root);

    const std::string verb = nuxsec::app::trim(args[0]);
    std::vector<std::string> rest;
    rest.reserve(args.size() > 0 ? args.size() - 1 : 0);
    for (size_t i = 1; i < args.size(); ++i)
    {
        rest.emplace_back(args[i]);
    }

    if (verb == "list")
    {
        if (!rest.empty())
        {
            throw std::runtime_error(kUsageMacro);
        }
        print_macro_list(std::cout, repo_root);
        return 0;
    }

    if (verb == "run")
    {
        if (rest.empty() || rest.size() > 2)
        {
            throw std::runtime_error(kUsageMacro);
        }
        const std::string macro_name = nuxsec::app::trim(rest[0]);
        const std::string call = (rest.size() == 2) ? nuxsec::app::trim(rest[1]) : "";

        const auto macro_path = resolve_macro_path(repo_root, macro_name);
        if (call.empty())
        {
            return run_root_macro_exec(repo_root, macro_path);
        }
        return run_root_macro_call(repo_root, macro_path, call);
    }

    if (rest.size() > 1)
    {
        throw std::runtime_error(kUsageMacro);
    }

    const std::string macro_name = verb;
    const std::string call = rest.empty() ? "" : nuxsec::app::trim(rest[0]);
    const auto macro_path = resolve_macro_path(repo_root, macro_name);
    if (call.empty())
    {
        return run_root_macro_exec(repo_root, macro_path);
    }
    return run_root_macro_call(repo_root, macro_path, call);
}

}

int main(int argc, char **argv)
{
    return nuxsec::app::run_guarded(
        [argc, argv]()
        {
            if (argc < 2)
            {
                print_main_help(std::cerr);
                return 1;
            }

            const std::string command = argv[1];
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv, 2);

            if (command == "help" || command == "-h" || command == "--help")
            {
                print_main_help(std::cout);
                return 0;
            }

            if (command == "art")
            {
                return run_art_command(args);
            }
            if (command == "sample")
            {
                return run_sample_command(args);
            }
            if (command == "event")
            {
                return run_event_command(args);
            }
            if (command == "macro")
            {
                return run_macro_command(args);
            }

            std::cerr << "Unknown command: " << command << "\n";
            print_main_help(std::cerr);
            return 1;
        });
}
