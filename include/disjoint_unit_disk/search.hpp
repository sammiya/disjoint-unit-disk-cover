#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "disjoint_unit_disk/int_scaled_point.hpp"

namespace disjoint_unit_disk {

enum class SearchStatus { kCoverable, kNonCoverable, kUndecided };

struct ScaledCenters {
  std::int64_t scale;
  std::vector<IntScaledPoint> centers;
};

struct SearchResult {
  SearchStatus status;
  std::optional<ScaledCenters>
      sample_centers; // only when status == kCoverable or kUndecided
};

SearchResult
Search(const std::vector<std::vector<IntScaledPoint>> &point_groups_scaled,
       std::int64_t scale);

} // namespace disjoint_unit_disk
