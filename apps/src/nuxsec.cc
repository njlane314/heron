/* -- C++ -- */
/**
 *  @file  apps/src/nuxsec.cc
 *
 *  @brief Unified CLI for Nuxsec utilities.
 */

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sys/wait.h>

#include <TROOT.h>
#include <TSystem.h>

#include "AppUtils.hh"
namespace
{

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
        const auto out = (repo_root / "build" / "plot").string();
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

int exec_root_macro(const std::filesystem::path &repo_root,
                    const std::filesystem::path &macro_path,
                    const std::string &call_cmd)
{
    ensure_plot_env(repo_root);
    add_plot_include_paths(repo_root);

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

int handle_macro_command(const std::vector<std::string> &args)
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
            return exec_root_macro(repo_root, macro_path, "");
        }
        return exec_root_macro(repo_root, macro_path, call);
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
        return exec_root_macro(repo_root, macro_path, "");
    }
    return exec_root_macro(repo_root, macro_path, call);
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

            const std::pair<const char *, const char *> driver_map[] = {
                {"art", "nuxsecArtFileIOdriver"},
                {"sample", "nuxsecSampleIOdriver"},
                {"event", "nuxsecEventIOdriver"}
            };
            for (const auto &entry : driver_map)
            {
                if (command == entry.first)
                {
                    return dispatch_driver_command(entry.second, args);
                }
            }
            if (command == "macro")
            {
                return handle_macro_command(args);
            }

            std::cerr << "Unknown command: " << command << "\n";
            print_main_help(std::cerr);
            return 1;
        });
}
