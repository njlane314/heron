/* -- C++ -- */
/**
 *  @file  plot/include/PlotEnv.hh
 *
 *  @brief Environment helpers for plot output locations and formats, handling
 *         repository-relative paths and default plotting conventions.
 *
 *  Conventions:
 *    - All *relative* plot outputs are rooted under HERON_PLOT_DIR.
 *    - Defaults:
 *        HERON_PLOT_DIR    = <repo>/scratch/plot
 *        HERON_PLOT_FORMAT = pdf
 *
 *  The CLI sets HERON_REPO_ROOT and (if unset) HERON_PLOT_DIR so plots land
 *  consistently regardless of where the binary is invoked from. When using
 *  workspace switching, the CLI derives HERON_PLOT_DIR from HERON_PLOT_BASE
 *  and HERON_SET.
 */

#ifndef HERON_PLOT_PLOT_ENV_H
#define HERON_PLOT_PLOT_ENV_H

#include <cstdlib>
#include <filesystem>
#include <string>
#include <string_view>
#include <system_error>


namespace nu
{

inline constexpr int kCanvasWidth = 1200;
inline constexpr int kCanvasHeight = 700;

inline const char *getenv_cstr(const char *name) noexcept
{
    return std::getenv(name);
}

inline std::filesystem::path repo_root_dir()
{
    if (const char *e = getenv_cstr("HERON_REPO_ROOT"))
    {
        return std::filesystem::path(e);
    }
    return std::filesystem::current_path();
}

inline std::filesystem::path plot_output_dir_path()
{
    if (const char *e = getenv_cstr("HERON_PLOT_DIR"))
    {
        return std::filesystem::path(e);
    }
    return repo_root_dir() / "scratch" / "plot";
}

inline std::string plot_output_dir()
{
    return plot_output_dir_path().string();
}

inline std::string plot_image_format()
{
    if (const char *e = getenv_cstr("HERON_PLOT_FORMAT"))
    {
        return std::string(e);
    }
    return "pdf";
}

inline std::filesystem::path release_dir_path()
{
    if (const char *e = getenv_cstr("HERON_RELEASE_DIR"))
    {
        return std::filesystem::path(e);
    }
    return repo_root_dir() / "release";
}

inline void ensure_directory(const std::filesystem::path &dir)
{
    if (dir.empty())
    {
        return;
    }
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
}

// Resolve an output file, *always* rooting relative paths under the plot dir.
// If no extension is provided, the plot format is used (default: pdf).
inline std::filesystem::path resolve_output_file(std::string_view user_path,
                                                 std::string_view default_base,
                                                 std::string_view default_ext = {})
{
    const std::string ext = default_ext.empty() ? plot_image_format() : std::string(default_ext);

    std::filesystem::path p;
    if (user_path.empty())
    {
        p = std::filesystem::path(default_base);
    }
    else
    {
        p = std::filesystem::path(std::string(user_path));
    }

    if (p.is_relative())
    {
        p = plot_output_dir_path() / p;
    }

    if (p.extension().empty())
    {
        p.replace_extension(ext);
    }

    ensure_directory(p.parent_path());
    return p;
}

inline std::filesystem::path plot_output_file(std::string_view base,
                                              std::string_view ext = {})
{
    return resolve_output_file({}, base, ext);
}

} // namespace nu


#endif
