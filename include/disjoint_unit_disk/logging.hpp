#pragma once

#include <string_view>

namespace disjoint_unit_disk::logging {

enum class Level {
  kInfo,
  kWarning,
};

void SetMinimumLevel(Level level);

void Log(Level level, std::string_view message);

} // namespace disjoint_unit_disk::logging
