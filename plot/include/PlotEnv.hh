/* -- C++ -- */
/**
 *  @file  plot/include/PlotEnv.hh
 *
 *  @brief Environment helpers for plot output locations and formats.
 *
 *  Conventions:
 *    - All *relative* plot outputs are rooted under NUXSEC_PLOT_DIR.
 *    - Defaults:
 *        NUXSEC_PLOT_DIR    = <repo>/build/plot
 *        NUXSEC_PLOT_FORMAT = pdf
 *
 *  The CLI sets NUXSEC_REPO_ROOT and (if unset) NUXSEC_PLOT_DIR so plots land
 *  consistently regardless of where the binary is invoked from.
 */

#ifndef Nuxsec_plot_PlotEnv_H_INCLUDED
#define Nuxsec_plot_PlotEnv_H_INCLUDED

#include <cstdlib>
#include <filesystem>
#include <string>
#include <string_view>
#include <system_error>

namespace nuxsec
{
namespace plot
{

inline constexpr int kCanvasWidth = 1200;
inline constexpr int kCanvasHeight = 700;

inline const char *getenv_cstr(const char *name) noexcept
{
    return std::getenv(name);
}

inline std::filesystem::path repo_root_dir()
{
    if (const char *e = getenv_cstr("NUXSEC_REPO_ROOT"))
    {
        return std::filesystem::path(e);
    }
    return std::filesystem::current_path();
}

inline std::filesystem::path plot_output_dir_path()
{
    if (const char *e = getenv_cstr("NUXSEC_PLOT_DIR"))
    {
        return std::filesystem::path(e);
    }
    return repo_root_dir() / "build" / "plot";
}

inline std::string plot_output_dir()
{
    return plot_output_dir_path().string();
}

inline std::string plot_image_format()
{
    if (const char *e = getenv_cstr("NUXSEC_PLOT_FORMAT"))
    {
        return std::string(e);
    }
    return "pdf";
}

inline std::filesystem::path release_dir_path()
{
    if (const char *e = getenv_cstr("NUXSEC_RELEASE_DIR"))
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

} // namespace plot
} // namespace nuxsec

#endif
