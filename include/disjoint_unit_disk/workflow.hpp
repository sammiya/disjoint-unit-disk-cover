#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "disjoint_unit_disk/int_scaled_point.hpp"
#include "disjoint_unit_disk/search.hpp"

namespace disjoint_unit_disk {

enum class CoverStatus { kCoverable, kNonCoverable, kUndecided };

struct DisjointCoverReport {
  CoverStatus status;
  std::int64_t scale;
  std::vector<IntScaledPoint> points_scaled;
  std::optional<size_t> partitions_evaluated;
  std::optional<std::vector<std::vector<size_t>>> sample_partition;
  std::optional<ScaledCenters> sample_centers;
};

DisjointCoverReport
AnalyzeDisjointUnitDiskCover(const std::vector<IntScaledPoint> &scaled_points,
                             std::int64_t scale);

DisjointCoverReport AnalyzeDisjointUnitDiskCoverWithGroups(
    const std::vector<IntScaledPoint> &scaled_points,
    const std::vector<size_t> &group_ids, std::int64_t scale);

} // namespace disjoint_unit_disk
