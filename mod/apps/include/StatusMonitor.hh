/* -- C++ -- */
#ifndef HERON_APPS_STATUS_MONITOR_H
#define HERON_APPS_STATUS_MONITOR_H

#include <chrono>
#include <condition_variable>
#include <ctime>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

#include "AppLog.hh"



class StatusMonitor
{
public:
    StatusMonitor(const std::string &log_prefix,
                  const std::string &message,
                  const std::chrono::seconds interval = std::chrono::minutes(1))
        : log_prefix_(log_prefix),
          message_(message),
          interval_(interval),
          start_time_(std::chrono::steady_clock::now()),
          worker_(&StatusMonitor::run_loop, this)
    {
    }

    ~StatusMonitor()
    {
        stop();
    }

    StatusMonitor(const StatusMonitor &) = delete;
    StatusMonitor &operator=(const StatusMonitor &) = delete;

    void stop()
    {
        {
            std::lock_guard<std::mutex> guard(mutex_);
            done_ = true;
        }
        cv_.notify_all();
        if (worker_.joinable())
        {
            worker_.join();
        }
    }

private:
    void run_loop()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        while (!done_)
        {
            if (cv_.wait_for(lock, interval_, [this]() { return done_; }))
            {
                break;
            }
            std::ostringstream out;
            out << message_
                << " time=" << format_timestamp()
                << " elapsed=" << format_elapsed_seconds() << "s"
                << " interval=" << interval_.count() << "s";
            log_info(log_prefix_, out.str());
        }
    }

    std::string format_timestamp() const
    {
        const auto now = std::chrono::system_clock::now();
        const std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm tm_snapshot;
        localtime_r(&now_time, &tm_snapshot);
        std::ostringstream out;
        out << std::put_time(&tm_snapshot, "%Y-%m-%dT%H:%M:%S%z");
        return out.str();
    }

    long long format_elapsed_seconds() const
    {
        const auto elapsed = std::chrono::steady_clock::now() - start_time_;
        return std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
    }

    std::string log_prefix_;
    std::string message_;
    std::chrono::seconds interval_;
    std::chrono::steady_clock::time_point start_time_;
    std::mutex mutex_;
    std::condition_variable cv_;
    bool done_ = false;
    std::thread worker_;
};



#endif // HERON_APPS_STATUS_MONITOR_H
