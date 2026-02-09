/* -- C++ -- */
/**
 *  @file  apps/src/nuxsec.cc
 *
 *  @brief Unified CLI for Nuxsec utilities.
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

#include <sys/wait.h>

#include <TInterpreter.h>
#include <TROOT.h>
#include <TSystem.h>

#include "AppUtils.hh"


const char *kUsageMacro =
    "Usage: nuxsec macro MACRO.C [CALL]\n"
    "       nuxsec macro list\n"
    "\nEnvironment:\n"
    "  NUXSEC_PLOT_BASE    Plot base directory (default: <repo>/scratch/plot)\n"
    "  NUXSEC_PLOT_DIR     Output directory override (default: NUXSEC_PLOT_BASE/<set>)\n"
    "  NUXSEC_PLOT_FORMAT  Output extension (default: pdf)\n"
    "  NUXSEC_SET          Workspace selector (default: template)\n";

const char *kMainBanner =
    "███╗   ██╗██╗   ██╗██╗  ██╗███████╗███████╗ ██████╗\n"
    "████╗  ██║██║   ██║╚██╗██╔╝██╔════╝██╔════╝██╔════╝\n"
    "██╔██╗ ██║██║   ██║ ╚███╔╝ ███████╗█████╗  ██║     \n"
    "██║╚██╗██║██║   ██║ ██╔██╗ ╚════██║██╔══╝  ██║     \n"
    "██║ ╚████║╚██████╔╝██╔╝ ██╗███████║███████╗╚██████╗\n"
    "╚═╝  ╚═══╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝╚══════╝ ╚═════╝\n";

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
    out << kMainBanner << "\n"
        << "Neutrino cross-section analysis CLI for provenance, samples, events,\n"
        << "and plots.\n\n"
        << "Usage: nuxsec <command> [args]\n\n"
        << "Commands:\n"
        << "  art         Aggregate art provenance for an input\n"
        << "  sample      Aggregate Sample ROOT files from art provenance\n"
        << "  event       Build event-level output from aggregated samples\n"
        << "  macro       Run ROOT macros (plotting or standalone)\n"
        << "  status      Log status for executable binaries\n"
        << "  paths       Print resolved workspace paths\n"
        << "  env         Print environment exports for a workspace\n"
        << "\nGlobal options:\n"
        << "  -S, --set   Workspace selector (default: template)\n"
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
    return path_from_env_or_default("NUXSEC_OUT_BASE",
                                    repo_root / "scratch" / "out");
}

std::filesystem::path plot_base_dir(const std::filesystem::path &repo_root)
{
    return path_from_env_or_default("NUXSEC_PLOT_BASE",
                                    repo_root / "scratch" / "plot");
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
    return path_from_env_or_default("NUXSEC_PLOT_DIR", out);
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

std::filesystem::path resolve_driver_path(const std::string &driver_name)
{
    std::vector<std::filesystem::path> candidates;

    if (const char *driver_dir = std::getenv("NUXSEC_DRIVER_DIR"))
    {
        candidates.emplace_back(driver_dir);
    }

    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec)
    {
        candidates.push_back(exe.parent_path());
    }

    const auto repo_root = find_repo_root();
    candidates.push_back(repo_root / "build" / "bin");

    for (const auto &base : candidates)
    {
        const auto candidate = base / driver_name;
        if (std::filesystem::exists(candidate))
        {
            return candidate;
        }
    }

    return driver_name;
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
    if (!gSystem->Getenv("NUXSEC_REPO_ROOT"))
    {
        gSystem->Setenv("NUXSEC_REPO_ROOT", repo_root.string().c_str());
    }
    if (!gSystem->Getenv("NUXSEC_PLOT_DIR"))
    {
        const auto out = plot_dir(repo_root).string();
        gSystem->Setenv("NUXSEC_PLOT_DIR", out.c_str());
    }
}

int dispatch_driver_command(const std::string &driver_name,
                            const std::vector<std::string> &args)
{
    const auto driver_path = resolve_driver_path(driver_name);
    if (std::filesystem::exists(driver_path) && !is_executable(driver_path))
    {
        throw std::runtime_error("Driver is not executable: " + driver_path.string());
    }
    std::ostringstream command;
    command << shell_quote(driver_path.string());
    for (const auto &arg : args)
    {
        command << " " << shell_quote(arg);
    }

    const int result = std::system(command.str().c_str());
    if (result == -1)
    {
        throw std::runtime_error("Failed to launch driver: " + driver_path.string());
    }

    if (WIFEXITED(result))
    {
        return WEXITSTATUS(result);
    }
    if (WIFSIGNALED(result))
    {
        return 128 + WTERMSIG(result);
    }

    return result;
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
        const auto standalone_candidate = repo_root / "standalone" / "macro" / candidate;
        if (std::filesystem::exists(standalone_candidate))
        {
            return standalone_candidate;
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
    add(repo_root / "apps" / "include");
}

void ensure_plot_lib_loaded(const std::filesystem::path &repo_root)
{
    const auto lib_dir = repo_root / "build" / "lib";
    if (std::filesystem::exists(lib_dir))
    {
        // Allow ROOT/cling to find project shared libraries when running macros.
        gSystem->AddDynamicPath(lib_dir.string().c_str());
    }

    // NOTE: nuxsec binary links IO + ANA, but not PLOT. Plot macros need this loaded.
    const auto plot_lib = lib_dir / "libNuxsecPlot.so";
    if (std::filesystem::exists(plot_lib))
    {
        const int rc = gSystem->Load(plot_lib.string().c_str());
        if (rc < 0)
        {
            throw std::runtime_error("Failed to load plot library: " + plot_lib.string());
        }
        return;
    }

    const int rc = gSystem->Load("libNuxsecPlot.so");
    if (rc < 0)
    {
        throw std::runtime_error("Failed to load plot library: libNuxsecPlot.so");
    }
}

void print_paths(std::ostream &out, const std::filesystem::path &repo_root)
{
    out << "NUXSEC_REPO_ROOT=" << repo_root.string() << "\n";
    out << "NUXSEC_SET=" << workspace_set() << "\n";
    out << "NUXSEC_OUT_BASE=" << out_base_dir(repo_root).string() << "\n";
    out << "NUXSEC_PLOT_BASE=" << plot_base_dir(repo_root).string() << "\n";
    out << "ART_DIR=" << stage_dir(repo_root, "NUXSEC_ART_DIR", "art").string() << "\n";
    out << "SAMPLE_DIR=" << stage_dir(repo_root, "NUXSEC_SAMPLE_DIR", "sample").string() << "\n";
    out << "EVENT_DIR=" << stage_dir(repo_root, "NUXSEC_EVENT_DIR", "event").string() << "\n";
    out << "PLOT_DIR=" << plot_dir(repo_root).string() << "\n";
}

int handle_paths_command(const std::vector<std::string> &args,
                         const std::filesystem::path &repo_root)
{
    if (!args.empty())
    {
        throw std::runtime_error("Usage: nuxsec paths");
    }
    print_paths(std::cout, repo_root);
    return 0;
}

int handle_env_command(const std::vector<std::string> &args,
                       const std::filesystem::path &repo_root)
{
    if (args.size() > 1)
    {
        throw std::runtime_error("Usage: nuxsec env [SET]");
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

    std::cout << "export NUXSEC_SET=" << shell_quote(set_value) << "\n";
    std::cout << "export NUXSEC_OUT_BASE=" << shell_quote(out_base_dir(repo_root).string()) << "\n";
    std::cout << "export NUXSEC_PLOT_BASE=" << shell_quote(plot_base_dir(repo_root).string()) << "\n";
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

    list_macros(repo_root / "plot" / "macro", "Plot macros in", "");
    list_macros(repo_root / "standalone" / "macro", "Standalone macros in", "");
    list_macros(repo_root / "evd" / "macro", "Event-display macros in", "evd/macro/");
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

        const auto macro_path = resolve_macro_path(repo_root, macro_name);
        if (call.empty())
        {
            return exec_root_macro(repo_root, macro_path, "");
        }
        return exec_root_macro(repo_root, macro_path, call);
    }

    if (rest.size() > 1)
    {
        throw std::runtime_error(kUsageMacro);
    }

    const std::string macro_name = verb;
    const std::string call = rest.empty() ? "" : trim(rest[0]);
    const auto macro_path = resolve_macro_path(repo_root, macro_name);
    if (call.empty())
    {
        return exec_root_macro(repo_root, macro_path, "");
    }
    return exec_root_macro(repo_root, macro_path, call);
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
        throw std::runtime_error("Usage: nuxsec status [--interval SECONDS] [--count COUNT] [--once]");
    }
    return opts;
}

std::vector<std::filesystem::path> status_dirs(const std::filesystem::path &repo_root)
{
    std::vector<std::filesystem::path> dirs;
    if (const char *driver_dir = std::getenv("NUXSEC_DRIVER_DIR"))
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
    log_info("nuxsec", start_message.str());

    long long completed = 0;
    while (opts.count == 0 || completed < opts.count)
    {
        const auto executables = collect_executables(repo_root);
        std::ostringstream summary;
        summary << "action=exe_status_scan status=complete executables="
                << format_count(
                       static_cast<long long>(executables.size()));
        log_info("nuxsec", summary.str());

        if (executables.empty())
        {
            log_warning(
                "nuxsec",
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
                log_info("nuxsec", message.str());
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
            std::cout << "Usage: nuxsec paths\n";
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
            std::cout << "Usage: nuxsec env [SET]\n";
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
            std::cout << "Usage: nuxsec status [--interval SECONDS] [--count COUNT] [--once]\n";
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
            return dispatch_driver_command("nuxsecArtFileIOdriver", args);
        },
        []()
        {
            std::cout << "Usage: nuxsec art <args>\n";
        }
    });
    table.push_back(CommandEntry{
        "sample",
        [](const std::vector<std::string> &args)
        {
            return dispatch_driver_command("nuxsecSampleIOdriver", args);
        },
        []()
        {
            std::cout << "Usage: nuxsec sample <args>\n";
        }
    });
    table.push_back(CommandEntry{
        "event",
        [repo_root](const std::vector<std::string> &args)
        {
            if (args.size() == 1 && !is_help_arg(args[0]))
            {
                std::vector<std::string> rewritten;
                rewritten.push_back(default_samples_tsv(repo_root).string());
                rewritten.push_back(args[0]);
                return dispatch_driver_command("nuxsecEventIOdriver", rewritten);
            }
            if (args.size() == 2 && !is_help_arg(args[0]) && !is_help_arg(args[1]))
            {
                if (has_suffix(args[0], ".root"))
                {
                    std::vector<std::string> rewritten;
                    rewritten.push_back(default_samples_tsv(repo_root).string());
                    rewritten.push_back(args[0]);
                    rewritten.push_back(args[1]);
                    return dispatch_driver_command("nuxsecEventIOdriver", rewritten);
                }
            }
            return dispatch_driver_command("nuxsecEventIOdriver", args);
        },
        []()
        {
            std::cout << "Usage: nuxsec event <args>\n";
        }
    });
    return table;
}

int main(int argc, char **argv)
{
    return run_guarded(
        "nuxsec",
        [argc, argv]()
        {
            int i = 1;
            const GlobalOptions global_opts = parse_global(i, argc, argv);
            if (!global_opts.set.empty())
            {
                ::setenv("NUXSEC_SET", global_opts.set.c_str(), 1);
            }

            if (i >= argc)
            {
                print_main_help(std::cerr);
                return 1;
            }

            const auto repo_root = find_repo_root();
            if (!getenv_cstr("NUXSEC_REPO_ROOT"))
            {
                ::setenv("NUXSEC_REPO_ROOT", repo_root.string().c_str(), 1);
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
