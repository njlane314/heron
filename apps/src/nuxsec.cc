/* -- C++ -- */
/**
 *  @file  apps/src/nuxsec.cc
 *
 *  @brief Unified CLI for Nuxsec utilities.
 */

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>

#include "Utils.hh"
#include "ArtCommand.hh"
#include "SampleCommand.hh"

namespace
{

const char *kUsageArt = "Usage: nuxsec a|art|artio|artio-aggregate NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]";
const char *kUsageSample = "Usage: nuxsec s|samp|sample|sample-aggregate NAME:FILELIST";
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
        << "  art         Aggregate art provenance for a stage\n"
        << "  samp        Aggregate Sample ROOT files from art provenance\n"
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

int run_artio_command(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageArt << "\n";
        return 0;
    }

    const nuxsec::app::ArtArgs art_args = nuxsec::app::parse_art_args(args, kUsageArt);
    return nuxsec::app::run_artio(art_args, "nuxsec artio-aggregate");
}

int run_sample_command(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageSample << "\n";
        return 0;
    }

    const nuxsec::app::SampleArgs sample_args = nuxsec::app::parse_sample_args(args, kUsageSample);
    return nuxsec::app::run_sample(sample_args, "nuxsec sample-aggregate");
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
        const std::string macro_spec = nuxsec::app::trim(rest[0]);
        const std::string call = (rest.size() == 2) ? nuxsec::app::trim(rest[1]) : "";

        const auto macro_path = resolve_macro_path(repo_root, macro_spec);
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

    const std::string macro_spec = verb;
    const std::string call = rest.empty() ? "" : nuxsec::app::trim(rest[0]);
    const auto macro_path = resolve_macro_path(repo_root, macro_spec);
    if (call.empty())
    {
        return run_root_macro_exec(repo_root, macro_path);
    }
    return run_root_macro_call(repo_root, macro_path, call);
}

}

int main(int argc, char **argv)
{
    return nuxsec::app::run_with_exceptions(
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

            if (command == "artio-aggregate" || command == "artio" || command == "art" || command == "a")
            {
                return run_artio_command(args);
            }
            if (command == "sample-aggregate" || command == "sample" || command == "samp" || command == "s")
            {
                return run_sample_command(args);
            }
            if (command == "macro" || command == "m")
            {
                return run_macro_command(args);
            }

            std::cerr << "Unknown command: " << command << "\n";
            print_main_help(std::cerr);
            return 1;
        });
}
