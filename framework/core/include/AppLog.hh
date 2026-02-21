/* -- C++ -- */
/**
 *  @file  framework/core/include/AppLog.hh
 *
 *  @brief Logging helpers that format and emit consistent status messages for
 *         application command-line workflows, covering informational, warning,
 *         and error reporting paths.
 */
#ifndef HERON_CORE_APPLOG_H
#define HERON_CORE_APPLOG_H

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>




enum class Level
{
    kInfo,
    kSuccess,
    kWarn,
    kError
};

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
    out << prefix << " " << level_label << " ";
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




#endif // HERON_CORE_APPLOG_H
