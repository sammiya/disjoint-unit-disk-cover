#include "disjoint_unit_disk/logging.hpp"

#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>

namespace disjoint_unit_disk::logging {
namespace {
std::mutex log_mutex;
std::atomic<Level> minimum_level{Level::kWarning};

int LevelRank(Level level) {
  switch (level) {
  case Level::kInfo:
    return 0;
  case Level::kWarning:
    return 1;
  }
  return 1;
}

const char *ToPrefix(Level level) {
  switch (level) {
  case Level::kInfo:
    return "[INFO] ";
  case Level::kWarning:
    return "[WARN] ";
  }
  return "";
}

std::ostream &StreamFor(Level level) {
  if (level == Level::kWarning) {
    return std::cerr;
  }
  return std::clog;
}

std::string Timestamp() {
  using Clock = std::chrono::system_clock;
  const auto now = Clock::now();
  const auto time = Clock::to_time_t(now);
  const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(
                          now.time_since_epoch()) %
                      std::chrono::microseconds(1000000);
  std::tm tm_snapshot;
#if defined(_WIN32) || defined(_WIN64)
  localtime_s(&tm_snapshot, &time);
#else
  localtime_r(&time, &tm_snapshot);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm_snapshot, "%Y-%m-%d %H:%M:%S") << "."
      << std::setfill('0') << std::setw(6) << micros.count();
  return oss.str();
}
} // namespace

void SetMinimumLevel(Level level) {
  minimum_level.store(level, std::memory_order_relaxed);
}

void Log(Level level, std::string_view message) {
  if (LevelRank(level) <
      LevelRank(minimum_level.load(std::memory_order_relaxed))) {
    return;
  }
  std::lock_guard<std::mutex> guard(log_mutex);
  auto &stream = StreamFor(level);
  stream << Timestamp() << " " << ToPrefix(level) << message << std::endl;
}

} // namespace disjoint_unit_disk::logging
