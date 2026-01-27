/* -- C++ -- */
#ifndef NUXSEC_APPS_APPLOG_H
#define NUXSEC_APPS_APPLOG_H

#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <unistd.h>

namespace nuxsec
{

namespace app
{

namespace log
{

enum class Level
{
    kInfo,
    kSuccess,
    kWarn,
    kError
};

inline bool use_colour()
{
    static int enabled = -1;
    if (enabled >= 0)
    {
        return enabled != 0;
    }
    if (std::getenv("NO_COLOR") != nullptr ||
        std::getenv("NUXSEC_NO_COLOUR") != nullptr)
    {
        enabled = 0;
        return false;
    }
    enabled = (::isatty(::fileno(stderr)) != 0) ? 1 : 0;
    return enabled != 0;
}

inline const char *level_name(const Level level)
{
    switch (level)
    {
    case Level::kSuccess:
        return "DONE";
    case Level::kWarn:
        return "WARN";
    case Level::kError:
        return "ERROR";
    case Level::kInfo:
    default:
        return "INFO";
    }
}

inline const char *level_colour(const Level level)
{
    switch (level)
    {
    case Level::kSuccess:
        return "\033[1;32m";
    case Level::kWarn:
        return "\033[1;33m";
    case Level::kError:
        return "\033[1;31m";
    case Level::kInfo:
    default:
        return "\033[1;36m";
    }
}

inline std::string decorate(const std::string &text, const char *colour)
{
    if (!use_colour())
    {
        return text;
    }
    std::ostringstream out;
    out << colour << text << "\033[0m";
    return out.str();
}

inline std::string format_count(const long long count)
{
    std::ostringstream out;
    if (count >= 1000000)
    {
        out << std::fixed << std::setprecision(1)
            << (static_cast<double>(count) / 1000000.0) << "M";
    }
    else if (count >= 1000)
    {
        out << std::fixed << std::setprecision(1)
            << (static_cast<double>(count) / 1000.0) << "k";
    }
    else
    {
        out << count;
    }
    return out.str();
}

inline void log_line(const std::string &log_prefix,
                     const Level level,
                     const std::string &message)
{
    std::ostringstream out;
    const std::string prefix = "[" + log_prefix + "]";
    const std::string level_label = level_name(level);
    if (use_colour())
    {
        out << decorate(prefix, "\033[1;34m") << " "
            << decorate(level_label, level_colour(level)) << " ";
    }
    else
    {
        out << prefix << " " << level_label << " ";
    }
    out << message;
    std::cerr << out.str() << "\n";
}

inline void log_info(const std::string &log_prefix, const std::string &message)
{
    log_line(log_prefix, Level::kInfo, message);
}

inline void log_success(const std::string &log_prefix, const std::string &message)
{
    log_line(log_prefix, Level::kSuccess, message);
}

inline void log_warning(const std::string &log_prefix, const std::string &message)
{
    log_line(log_prefix, Level::kWarn, message);
}

inline void log_error(const std::string &log_prefix, const std::string &message)
{
    log_line(log_prefix, Level::kError, message);
}

inline void log_stage(const std::string &log_prefix,
                      const std::string &stage,
                      const std::string &detail = "")
{
    std::ostringstream out;
    out << "stage=" << stage;
    if (!detail.empty())
    {
        out << " " << detail;
    }
    log_info(log_prefix, out.str());
}

} // namespace log

} // namespace app

} // namespace nuxsec

#endif // NUXSEC_APPS_APPLOG_H
