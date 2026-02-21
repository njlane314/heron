/* -- C++ -- */
/**
 *  @file  framework/core/src/heron.cc
 *
 *  @brief Unified CLI for HERON utilities.
 */

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <TInterpreter.h>
#include <TROOT.h>
#include <TSystem.h>

#include "ArtCLI.hh"
#include "EventCLI.hh"
#include "AppUtils.hh"
#include "SampleCLI.hh"


const char *kUsageMacro =
    "Usage: heron macro MACRO.C [CALL]\n"
    "       heron macro list\n"
    "\nEnvironment:\n"
    "  HERON_MACRO_LIBRARY_DIR  In-repo macro library directory (default: <repo>/macros/library)\n"
    "  HERON_MACRO_PATH         Colon-separated extra macro directories (searched after library)\n"
    "  HERON_PLOT_BASE    Plot base directory (default: <repo>/scratch/plot)\n"
    "  HERON_PLOT_DIR     Output directory override (default: HERON_PLOT_BASE/<set>)\n"
    "  HERON_PLOT_FORMAT  Output extension (default: pdf)\n"
    "  HERON_SET          Workspace selector (default: out)\n";

bool is_help_arg(const std::string &arg)
{
    return arg == "-h" || arg == "--help";
}

bool has_suffix(const std::string &value, const std::string &suffix)
{
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

struct GlobalOptions
{
    std::string set;
};

struct CommandEntry
{
    const char *name;
    std::function<int(const std::vector<std::string> &)> handler;
    std::function<void()> help;
};

GlobalOptions parse_global(int &i, int argc, char **argv)
{
    GlobalOptions opt;
    for (; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--set" || arg == "-S")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("Missing value for --set");
            }
            opt.set = argv[++i];
            continue;
        }
        if (arg == "--")
        {
            ++i;
            break;
        }
        if (!arg.empty() && arg[0] == '-')
        {
            break;
        }
        break;
    }
    return opt;
}

void print_main_help(std::ostream &out)
{
    out << "HERON — Histogram and Event Relay for Orchestrated Normalisation.\n\n"
        << "Usage: heron <command> [args]\n\n"
        << "Commands:\n"
        << "  art         Aggregate art provenance for an input\n"
        << "  sample      Aggregate Sample ROOT files from art provenance\n"
        << "  event       Build event-level output from aggregated samples\n"
        << "  macro       Run ROOT macros (plotting or standalone)\n"
        << "  status      Log status for executable binaries\n"
        << "  paths       Print resolved workspace paths\n"
        << "  env         Print environment exports for a workspace\n"
        << "\nGlobal options:\n"
        << "  -S, --set   Workspace selector (default: out)\n"
        << "\nRun 'heron <command> --help' for command-specific usage.\n";
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
            if (std::filesystem::exists(base / "macros/.plot_driver.retired") ||
                std::filesystem::exists(base / "framework/modules/plot/macro/.plot_driver.retired"))
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

std::string string_from_env_or_default(const char *name,
                                       const std::string &fallback)
{
    if (const char *value = getenv_cstr(name))
    {
        return std::string(value);
    }
    return fallback;
}

std::filesystem::path path_from_env_or_default(
    const char *name,
    const std::filesystem::path &fallback)
{
    if (const char *value = getenv_cstr(name))
    {
        return std::filesystem::path(value);
    }
    return fallback;
}

std::filesystem::path out_base_dir(const std::filesystem::path &repo_root)
{
    const auto fallback = path_from_env_or_default("HERON_OUTPUT_DIR",
                                                   repo_root / "scratch" / "out");
    return path_from_env_or_default("HERON_OUT_BASE", fallback);
}

std::filesystem::path plot_base_dir(const std::filesystem::path &repo_root)
{
    return path_from_env_or_default("HERON_PLOT_BASE",
                                    repo_root / "scratch" / "plot");
}

std::filesystem::path macro_repo_dir(const std::filesystem::path &repo_root)
{
    const auto macro_submodule = repo_root / "macros";
    if (std::filesystem::exists(macro_submodule))
    {
        return macro_submodule;
    }
    return repo_root;
}

std::filesystem::path stage_dir(const std::filesystem::path &repo_root,
                                const char *override_env,
                                const std::string &stage)
{
    const auto fallback = out_base_dir(repo_root) / workspace_set() / stage;
    return path_from_env_or_default(override_env, fallback);
}

std::filesystem::path plot_dir(const std::filesystem::path &repo_root)
{
    std::filesystem::path out = plot_base_dir(repo_root);
    const std::string set = workspace_set();
    if (!set.empty())
    {
        out /= set;
    }
    return path_from_env_or_default("HERON_PLOT_DIR", out);
}

std::filesystem::path default_samples_tsv(const std::filesystem::path &repo_root)
{
    return out_base_dir(repo_root) / workspace_set() / "sample" / "samples.tsv";
}

std::string shell_quote(const std::string &value)
{
    if (value.empty())
    {
        return "''";
    }
    std::string quoted;
    quoted.reserve(value.size() + 2);
    quoted.push_back('\'');
    for (char c : value)
    {
        if (c == '\'')
        {
            quoted.append("'\\''");
        }
        else
        {
            quoted.push_back(c);
        }
    }
    quoted.push_back('\'');
    return quoted;
}

bool is_executable(const std::filesystem::path &path)
{
    std::error_code ec;
    const auto status = std::filesystem::status(path, ec);
    if (ec || !std::filesystem::is_regular_file(status))
    {
        return false;
    }
    const auto exec_perms = std::filesystem::perms::owner_exec |
                            std::filesystem::perms::group_exec |
                            std::filesystem::perms::others_exec;
    return (status.permissions() & exec_perms) != std::filesystem::perms::none;
}

void ensure_plot_env(const std::filesystem::path &repo_root)
{
    if (!gSystem->Getenv("HERON_REPO_ROOT"))
    {
        gSystem->Setenv("HERON_REPO_ROOT", repo_root.string().c_str());
    }
    if (!gSystem->Getenv("HERON_PLOT_DIR"))
    {
        const auto out = plot_dir(repo_root).string();
        gSystem->Setenv("HERON_PLOT_DIR", out.c_str());
    }
}

std::vector<std::string> split_path_list(const std::string &raw)
{
    std::vector<std::string> entries;
    std::stringstream ss(raw);
    std::string item;
    while (std::getline(ss, item, ':'))
    {
        const std::string value = trim(item);
        if (!value.empty())
        {
            entries.push_back(value);
        }
    }
    return entries;
}

std::filesystem::path macro_library_dir(const std::filesystem::path &repo_root)
{
    const std::filesystem::path fallback = macro_repo_dir(repo_root) / "library";
    return path_from_env_or_default("HERON_MACRO_LIBRARY_DIR", fallback);
}

std::vector<std::filesystem::path> macro_search_dirs(const std::filesystem::path &repo_root)
{
    auto append_if_exists = [](std::vector<std::filesystem::path> &dirs,
                               const std::filesystem::path &path)
    {
        if (!path.empty() && std::filesystem::exists(path))
        {
            dirs.push_back(path);
        }
    };

    std::vector<std::filesystem::path> dirs;

    append_if_exists(dirs, macro_library_dir(repo_root));

    const std::string raw_dirs = string_from_env_or_default("HERON_MACRO_PATH", "");
    const auto dir_entries = split_path_list(raw_dirs);
    for (const auto &entry : dir_entries)
    {
        std::filesystem::path p(entry);
        if (p.is_relative())
        {
            p = repo_root / p;
        }
        append_if_exists(dirs, p);
    }

    return dirs;
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

        const auto external_macro_dirs = macro_search_dirs(repo_root);
        for (const auto &macro_dir : external_macro_dirs)
        {
            const auto custom_candidate = macro_dir / candidate;
            if (std::filesystem::exists(custom_candidate))
            {
                return custom_candidate;
            }
        }

        const auto macro_root = macro_repo_dir(repo_root);
        const auto standalone_candidate = macro_root / "standalone" / candidate;
        if (std::filesystem::exists(standalone_candidate))
        {
            return standalone_candidate;
        }
        const auto macro_candidate = macro_root / candidate;
        if (std::filesystem::exists(macro_candidate))
        {
            return macro_candidate;
        }
    }
    return candidate;
}

void add_plot_include_paths(const std::filesystem::path &repo_root)
{
    auto add = [&](const std::filesystem::path &p)
    {
        const std::string include_flag = "-I" + p.string();
        gSystem->AddIncludePath(include_flag.c_str());
        if (gInterpreter)
        {
            gInterpreter->AddIncludePath(p.string().c_str());
        }
    };
    add(repo_root / "plot" / "include");
    add(repo_root / "ana" / "include");
    add(repo_root / "io" / "include");
    add(repo_root / "core" / "include");
}

void ensure_plot_lib_loaded(const std::filesystem::path &repo_root)
{
    const auto lib_dir = repo_root / "build" / "lib";
    if (std::filesystem::exists(lib_dir))
    {
        // Allow ROOT/cling to find project shared libraries when running macros.
        gSystem->AddDynamicPath(lib_dir.string().c_str());
    }

    // NOTE: heron binary links IO + ANA, but not PLOT. Plot macros need this loaded.
    const auto plot_lib = lib_dir / "libheronPlot.so";
    if (std::filesystem::exists(plot_lib))
    {
        const int rc = gSystem->Load(plot_lib.string().c_str());
        if (rc < 0)
        {
            throw std::runtime_error("Failed to load plot library: " + plot_lib.string());
        }
        return;
    }

    const int rc = gSystem->Load("libheronPlot.so");
    if (rc < 0)
    {
        throw std::runtime_error("Failed to load plot library: libheronPlot.so");
    }
}

void print_paths(std::ostream &out, const std::filesystem::path &repo_root)
{
    out << "HERON_REPO_ROOT=" << repo_root.string() << "\n";
    out << "HERON_SET=" << workspace_set() << "\n";
    out << "HERON_OUT_BASE=" << out_base_dir(repo_root).string() << "\n";
    out << "HERON_PLOT_BASE=" << plot_base_dir(repo_root).string() << "\n";
    out << "ART_DIR=" << stage_dir(repo_root, "HERON_ART_DIR", "art").string() << "\n";
    out << "SAMPLE_DIR=" << stage_dir(repo_root, "HERON_SAMPLE_DIR", "sample").string() << "\n";
    out << "EVENT_DIR=" << stage_dir(repo_root, "HERON_EVENT_DIR", "event").string() << "\n";
    out << "PLOT_DIR=" << plot_dir(repo_root).string() << "\n";
}

int handle_paths_command(const std::vector<std::string> &args,
                         const std::filesystem::path &repo_root)
{
    if (!args.empty())
    {
        throw std::runtime_error("Usage: heron paths");
    }
    print_paths(std::cout, repo_root);
    return 0;
}

int handle_env_command(const std::vector<std::string> &args,
                       const std::filesystem::path &repo_root)
{
    if (args.size() > 1)
    {
        throw std::runtime_error("Usage: heron env [SET]");
    }

    std::string set_value = workspace_set();
    if (!args.empty())
    {
        set_value = trim(args[0]);
    }

    if (set_value.empty())
    {
        throw std::runtime_error("Missing workspace set value");
    }

    std::cout << "export HERON_SET=" << shell_quote(set_value) << "\n";
    std::cout << "export HERON_OUT_BASE=" << shell_quote(out_base_dir(repo_root).string()) << "\n";
    std::cout << "export HERON_PLOT_BASE=" << shell_quote(plot_base_dir(repo_root).string()) << "\n";
    return 0;
}

bool is_plot_macro(const std::filesystem::path &repo_root,
                   const std::filesystem::path &macro_path)
{
    std::error_code ec;
    const auto rel = std::filesystem::relative(macro_path, repo_root, ec);
    if (ec)
    {
        return false;
    }
    if (rel.empty())
    {
        return false;
    }
    auto it = rel.begin();
    if (it == rel.end() || *it != "plot")
    {
        return false;
    }
    ++it;
    return it != rel.end() && *it == "macro";
}

int exec_root_macro(const std::filesystem::path &repo_root,
                    const std::filesystem::path &macro_path,
                    const std::string &call_cmd)
{
    ensure_plot_env(repo_root);
    add_plot_include_paths(repo_root);
    if (is_plot_macro(repo_root, macro_path))
    {
        ensure_plot_lib_loaded(repo_root);
    }

    if (!std::filesystem::exists(macro_path))
    {
        throw std::runtime_error("Macro not found at " + macro_path.string());
    }

    const bool has_call = !call_cmd.empty();
    if (has_call)
    {
        const std::string load_cmd = ".L " + macro_path.string();
        gROOT->ProcessLine(load_cmd.c_str());
        const long result = gROOT->ProcessLine(call_cmd.c_str());
        return static_cast<int>(result);
    }
    const std::string exec_cmd = ".x " + macro_path.string();
    const long result = gROOT->ProcessLine(exec_cmd.c_str());
    return static_cast<int>(result);
}

void print_macro_list(std::ostream &out, const std::filesystem::path &repo_root)
{
    auto list_macros = [&](const std::filesystem::path &dir,
                           const std::string &label,
                           const std::string &prefix)
    {
        out << label << " " << dir.string() << ":\n";
        if (!std::filesystem::exists(dir))
        {
            out << "  (none; directory not found)\n";
            return;
        }

        std::vector<std::string> macros;
        for (const auto &entry : std::filesystem::directory_iterator(dir))
        {
            if (!entry.is_regular_file())
            {
                continue;
            }
            const auto &path = entry.path();
            if (path.extension() == ".C")
            {
                macros.push_back(prefix + path.filename().string());
            }
        }

        std::sort(macros.begin(), macros.end());
        for (const auto &macro : macros)
        {
            out << "  " << macro << "\n";
        }
    };

    const auto external_dirs = macro_search_dirs(repo_root);
    for (const auto &dir : external_dirs)
    {
        list_macros(dir, "External macros in", "");
    }

    const auto macro_root = macro_repo_dir(repo_root);
    list_macros(macro_root, "Macros in", "");
    list_macros(macro_root / "standalone", "Standalone macros in", "standalone/");
}

int handle_macro_command(const std::vector<std::string> &args)
{
    if (args.empty())
    {
        std::cout << kUsageMacro << "\n";
        print_macro_list(std::cout, find_repo_root());
        return 0;
    }

    const auto repo_root = find_repo_root();
    ensure_plot_env(repo_root);

    const std::string verb = trim(args[0]);
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
        const std::string macro_name = trim(rest[0]);
        const std::string call = (rest.size() == 2) ? trim(rest[1]) : "";

        const std::filesystem::path macro_path = resolve_macro_path(repo_root, macro_name);
        return exec_root_macro(repo_root, macro_path, call);
    }


    if (rest.size() > 1)
    {
        throw std::runtime_error(kUsageMacro);
    }

    const std::string macro_name = verb;
    std::string call = rest.empty() ? "" : trim(rest[0]);
    const std::filesystem::path macro_path = resolve_macro_path(repo_root, macro_name);
    return exec_root_macro(repo_root, macro_path, call);
}


int handle_art_command(const std::vector<std::string> &args)
{
    return run_guarded(
        "heronArtFileIOdriver",
        [&]()
        {
            const ArtArgs art_args =
                parse_art_args(
                    args,
                    "Usage: heron art INPUT_NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]");
            return run(art_args, "heronArtFileIOdriver");
        });
}

int handle_sample_command(const std::vector<std::string> &args)
{
    return run_guarded(
        "heronSampleIOdriver",
        [&]()
        {
            const SampleArgs sample_args =
                parse_sample_args(
                    args,
                    "Usage: heron sample NAME:FILELIST");
            return run(sample_args, "heronSampleIOdriver");
        });
}

int handle_event_command(const std::vector<std::string> &args,
                         const std::filesystem::path &repo_root)
{
    return run_guarded(
        "heronEventIOdriver",
        [&]()
        {
            std::vector<std::string> rewritten = args;
            if (args.size() == 1)
            {
                rewritten.clear();
                rewritten.push_back(default_samples_tsv(repo_root).string());
                rewritten.push_back(args[0]);
            }
            else if (args.size() == 2 && has_suffix(args[0], ".root"))
            {
                rewritten.clear();
                rewritten.push_back(default_samples_tsv(repo_root).string());
                rewritten.push_back(args[0]);
                rewritten.push_back(args[1]);
            }

            const EventArgs event_args =
                parse_event_args(
                    rewritten,
                    "Usage: heron event SAMPLE_LIST.tsv OUTPUT.root [SELECTION] [COLUMNS.tsv]");
            return run(event_args, "heronEventIOdriver");
        });
}

struct StatusOptions
{
    int interval_seconds = 60;
    long long count = 0;
};

StatusOptions parse_status_args(const std::vector<std::string> &args)
{
    StatusOptions opts;
    for (size_t i = 0; i < args.size(); ++i)
    {
        const std::string arg = trim(args[i]);
        if (arg == "--interval" || arg == "-i")
        {
            if (i + 1 >= args.size())
            {
                throw std::runtime_error("Missing value for --interval");
            }
            opts.interval_seconds = std::stoi(args[++i]);
            if (opts.interval_seconds <= 0)
            {
                throw std::runtime_error("Interval must be positive");
            }
            continue;
        }
        if (arg == "--count" || arg == "-n")
        {
            if (i + 1 >= args.size())
            {
                throw std::runtime_error("Missing value for --count");
            }
            opts.count = std::stoll(args[++i]);
            if (opts.count <= 0)
            {
                throw std::runtime_error("Count must be positive");
            }
            continue;
        }
        if (arg == "--once")
        {
            opts.count = 1;
            continue;
        }
        throw std::runtime_error("Usage: heron status [--interval SECONDS] [--count COUNT] [--once]");
    }
    return opts;
}

std::vector<std::filesystem::path> status_dirs(const std::filesystem::path &repo_root)
{
    std::vector<std::filesystem::path> dirs;
    if (const char *driver_dir = std::getenv("HERON_DRIVER_DIR"))
    {
        dirs.emplace_back(driver_dir);
    }
    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec)
    {
        dirs.push_back(exe.parent_path());
    }
    dirs.push_back(repo_root / "build" / "bin");

    std::vector<std::filesystem::path> unique_dirs;
    std::set<std::filesystem::path> seen;
    for (const auto &dir : dirs)
    {
        if (dir.empty())
        {
            continue;
        }
        if (!seen.insert(dir).second)
        {
            continue;
        }
        unique_dirs.push_back(dir);
    }
    return unique_dirs;
}

std::vector<std::filesystem::path> collect_executables(
    const std::filesystem::path &repo_root)
{
    std::vector<std::filesystem::path> executables;
    const auto dirs = status_dirs(repo_root);
    for (const auto &dir : dirs)
    {
        if (!std::filesystem::exists(dir))
        {
            continue;
        }
        if (!std::filesystem::is_directory(dir))
        {
            if (is_executable(dir))
            {
                executables.push_back(dir);
            }
            continue;
        }
        for (const auto &entry : std::filesystem::directory_iterator(dir))
        {
            if (!entry.is_regular_file())
            {
                continue;
            }
            const auto &path = entry.path();
            if (is_executable(path))
            {
                executables.push_back(path);
            }
        }
    }
    std::sort(executables.begin(), executables.end());
    executables.erase(std::unique(executables.begin(), executables.end()),
                      executables.end());
    return executables;
}

int handle_status_command(const std::vector<std::string> &args,
                          const std::filesystem::path &repo_root)
{
    const StatusOptions opts = parse_status_args(args);
    std::ostringstream start_message;
    start_message << "action=exe_status_monitor status=start interval="
                  << opts.interval_seconds << "s";
    if (opts.count > 0)
    {
        start_message << " count=" << opts.count;
    }
    log_info("heron", start_message.str());

    long long completed = 0;
    while (opts.count == 0 || completed < opts.count)
    {
        const auto executables = collect_executables(repo_root);
        std::ostringstream summary;
        summary << "action=exe_status_scan status=complete executables="
                << format_count(
                       static_cast<long long>(executables.size()));
        log_info("heron", summary.str());

        if (executables.empty())
        {
            log_warning(
                "heron",
                "action=exe_status status=empty message=No executables found");
        }
        else
        {
            for (const auto &path : executables)
            {
                std::ostringstream message;
                message << "action=exe_status status=ok exe="
                        << path.filename().string()
                        << " path=" << path.string();
                log_info("heron", message.str());
            }
        }

        ++completed;
        if (opts.count != 0 && completed >= opts.count)
        {
            break;
        }
        std::this_thread::sleep_for(
            std::chrono::seconds(opts.interval_seconds));
    }
    return 0;
}

std::vector<CommandEntry> build_command_table(const std::filesystem::path &repo_root)
{
    std::vector<CommandEntry> table;
    table.push_back(CommandEntry{
        "help",
        [](const std::vector<std::string> &)
        {
            print_main_help(std::cout);
            return 0;
        },
        []()
        {
            print_main_help(std::cout);
        }
    });
    table.push_back(CommandEntry{
        "-h",
        [](const std::vector<std::string> &)
        {
            print_main_help(std::cout);
            return 0;
        },
        []()
        {
            print_main_help(std::cout);
        }
    });
    table.push_back(CommandEntry{
        "--help",
        [](const std::vector<std::string> &)
        {
            print_main_help(std::cout);
            return 0;
        },
        []()
        {
            print_main_help(std::cout);
        }
    });
    table.push_back(CommandEntry{
        "paths",
        [repo_root](const std::vector<std::string> &args)
        {
            return handle_paths_command(args, repo_root);
        },
        []()
        {
            std::cout << "Usage: heron paths\n";
        }
    });
    table.push_back(CommandEntry{
        "env",
        [repo_root](const std::vector<std::string> &args)
        {
            return handle_env_command(args, repo_root);
        },
        []()
        {
            std::cout << "Usage: heron env [SET]\n";
        }
    });
    table.push_back(CommandEntry{
        "status",
        [repo_root](const std::vector<std::string> &args)
        {
            return handle_status_command(args, repo_root);
        },
        []()
        {
            std::cout << "Usage: heron status [--interval SECONDS] [--count COUNT] [--once]\n";
        }
    });
    table.push_back(CommandEntry{
        "macro",
        [](const std::vector<std::string> &args)
        {
            return handle_macro_command(args);
        },
        []()
        {
            std::cout << kUsageMacro << "\n";
            print_macro_list(std::cout, find_repo_root());
        }
    });
    table.push_back(CommandEntry{
        "art",
        [](const std::vector<std::string> &args)
        {
            return handle_art_command(args);
        },
        []()
        {
            std::cout << "Usage: heron art INPUT_NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]\n";
        }
    });
    table.push_back(CommandEntry{
        "sample",
        [](const std::vector<std::string> &args)
        {
            return handle_sample_command(args);
        },
        []()
        {
            std::cout << "Usage: heron sample NAME:FILELIST\n";
        }
    });
    table.push_back(CommandEntry{
        "event",
        [repo_root](const std::vector<std::string> &args)
        {
            return handle_event_command(args, repo_root);
        },
        []()
        {
            std::cout << "Usage: heron event SAMPLE_LIST.tsv OUTPUT.root [SELECTION] [COLUMNS.tsv]\n";
        }
    });
    return table;
}

int main(int argc, char **argv)
{
    return run_guarded(
        "heron",
        [argc, argv]()
        {
            int i = 1;
            const GlobalOptions global_opts = parse_global(i, argc, argv);
            if (!global_opts.set.empty())
            {
                ::setenv("HERON_SET", global_opts.set.c_str(), 1);
            }

            if (i >= argc)
            {
                print_main_help(std::cerr);
                return 1;
            }

            const auto repo_root = find_repo_root();
            if (!getenv_cstr("HERON_REPO_ROOT"))
            {
                ::setenv("HERON_REPO_ROOT", repo_root.string().c_str(), 1);
            }

            const std::string command = argv[i++];
            const std::vector<std::string> args = collect_args(argc, argv, i);

            const auto command_table = build_command_table(repo_root);
            for (const auto &entry : command_table)
            {
                if (command == entry.name)
                {
                    if (!args.empty() && is_help_arg(args[0]))
                    {
                        entry.help();
                        return 0;
                    }
                    return entry.handler(args);
                }
            }

            std::cerr << "Unknown command: " << command << "\n";
            print_main_help(std::cerr);
            return 1;
        });
}
