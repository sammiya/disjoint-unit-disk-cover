#include "disjoint_unit_disk/workflow.hpp"
#include "disjoint_unit_disk/enumeration.hpp"

#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "disjoint_unit_disk/logging.hpp"
namespace disjoint_unit_disk {

namespace {

bool WillAddOverflow(std::int64_t a, std::int64_t b) {
  if (b > 0 && a > std::numeric_limits<std::int64_t>::max() - b) {
    return true;
  }
  if (b < 0 && a < std::numeric_limits<std::int64_t>::min() - b) {
    return true;
  }
  return false;
}

bool WillSubOverflow(std::int64_t a, std::int64_t b) {
  // Checks whether (a - b) would overflow int64_t without computing -b.
  if (b > 0 && a < std::numeric_limits<std::int64_t>::min() + b) {
    return true;
  }
  if (b < 0 && a > std::numeric_limits<std::int64_t>::max() + b) {
    return true;
  }
  return false;
}

bool WillSquareOverflow(std::int64_t v) {
  if (v == 0) {
    return false;
  }
  if (v == std::numeric_limits<std::int64_t>::min()) {
    return true;
  }
  const std::int64_t abs_v = v < 0 ? -v : v;
  return abs_v > std::numeric_limits<std::int64_t>::max() / abs_v;
}

std::int64_t SafeAdd(std::int64_t a, std::int64_t b) {
  if (WillAddOverflow(a, b)) {
    throw std::runtime_error(
        "Input coordinates overflow when expanded by scale");
  }
  return a + b;
}

std::int64_t SafeSub(std::int64_t a, std::int64_t b) {
  if (WillSubOverflow(a, b)) {
    throw std::runtime_error(
        "Input coordinates overflow when expanded by scale");
  }
  return a - b;
}

void ValidateOverflowFreeExpansion(const std::vector<IntScaledPoint> &points,
                                   std::int64_t scale) {
  if (scale < 1) {
    throw std::invalid_argument("scale must be an integer >= 1");
  }

  // Verify scale^2 and 4*scale^2 fit in int64_t before multiplying.
  const std::int64_t max = std::numeric_limits<std::int64_t>::max();
  if (scale > max / scale) {
    throw std::runtime_error("scale is too large: scale^2 overflows int64_t");
  }
  const std::int64_t scale_sq = scale * scale; // safe after check above
  if (scale_sq > max / 4) {
    throw std::runtime_error(
        "scale is too large: 4*scale^2 overflows int64_t");
  }

  // Each point expands to a rectangle with corners at (x±scale, y±scale); later
  // stages compute squared distances between such corners. If any corner
  // coordinate overflows on expansion, or if any squared distance between any
  // pair of corners would overflow, bail out early.
  std::vector<IntScaledPoint> corners;
  corners.reserve(points.size() * 4);

  for (const auto &p : points) {
    const std::int64_t x_plus = SafeAdd(p.x, scale);
    const std::int64_t x_minus = SafeSub(p.x, scale);
    const std::int64_t y_plus = SafeAdd(p.y, scale);
    const std::int64_t y_minus = SafeSub(p.y, scale);
    corners.push_back(IntScaledPoint{x_minus, y_minus});
    corners.push_back(IntScaledPoint{x_minus, y_plus});
    corners.push_back(IntScaledPoint{x_plus, y_minus});
    corners.push_back(IntScaledPoint{x_plus, y_plus});
  }

  for (size_t i = 0; i < corners.size(); ++i) {
    for (size_t j = i + 1; j < corners.size(); ++j) {
      const std::int64_t dx = SafeSub(corners[i].x, corners[j].x);
      const std::int64_t dy = SafeSub(corners[i].y, corners[j].y);
      if (WillSquareOverflow(dx) || WillSquareOverflow(dy)) {
        throw std::runtime_error(
            "Input coordinates overflow when computing squared distances");
      }
      const std::int64_t dx_sq = dx * dx;
      const std::int64_t dy_sq = dy * dy;
      if (WillAddOverflow(dx_sq, dy_sq)) {
        throw std::runtime_error(
            "Input coordinates overflow when computing squared distances");
      }
    }
  }
}

template <typename T> void EnsureNonEmptyPoints(const std::vector<T> &points) {
  if (points.empty()) {
    throw std::invalid_argument(
        "at least one point is required to analyze disjoint covers");
  }
}

std::vector<std::vector<IntScaledPoint>>
BuildGroups(const std::vector<IntScaledPoint> &points,
            const std::vector<std::vector<size_t>> &partition) {
  std::vector<std::vector<IntScaledPoint>> groups;
  groups.reserve(partition.size());
  for (const auto &indices : partition) {
    std::vector<IntScaledPoint> group;
    group.reserve(indices.size());
    for (size_t index : indices) {
      if (index >= points.size()) {
        throw std::out_of_range("partition index out of range");
      }
      group.push_back(points[index]);
    }
    groups.push_back(std::move(group));
  }
  return groups;
}

SearchResult SearchPartition(const std::vector<IntScaledPoint> &points,
                             const std::vector<std::vector<size_t>> &partition,
                             std::int64_t scale) {
  const auto groups = BuildGroups(points, partition);
  return Search(groups, scale);
}

} // namespace

DisjointCoverReport
AnalyzeDisjointUnitDiskCover(const std::vector<IntScaledPoint> &scaled_points,
                             std::int64_t scale) {
  EnsureNonEmptyPoints(scaled_points);
  ValidateOverflowFreeExpansion(scaled_points, scale);

  const auto cover_sets = EnumerateDiskCovers(scaled_points, scale);
  auto partitions = EnumerateCoverPartitions(cover_sets, scaled_points.size());

  logging::Log(logging::Level::kInfo, "Evaluating partition search across " +
                                          std::to_string(partitions.size()) +
                                          " candidates");

  std::optional<std::vector<std::vector<size_t>>> sample_partition;
  std::optional<ScaledCenters> sample_centers;
  CoverStatus status = CoverStatus::kNonCoverable;

  for (size_t partition_idx = 0; partition_idx < partitions.size();
       ++partition_idx) {
    logging::Log(logging::Level::kInfo,
                 "Analyzing partition " + std::to_string(partition_idx + 1) +
                     "/" + std::to_string(partitions.size()));
    const auto &partition = partitions[partition_idx];
    const auto result = SearchPartition(scaled_points, partition, scale);
    if (result.status == SearchStatus::kCoverable) {
      sample_partition = partition;
      sample_centers = result.sample_centers;
      status = CoverStatus::kCoverable;
      logging::Log(logging::Level::kInfo,
                   "Found disjoint cover on partition " +
                       std::to_string(partition_idx + 1));
      break;
    } else if (result.status == SearchStatus::kUndecided) {
      status = CoverStatus::kUndecided;
    }
  }

  return DisjointCoverReport{
      .status = status,
      .scale = scale,
      .points_scaled = scaled_points,
      .partitions_evaluated = partitions.size(),
      .sample_partition = sample_partition,
      .sample_centers = sample_centers,
  };
}

DisjointCoverReport AnalyzeDisjointUnitDiskCoverWithGroups(
    const std::vector<IntScaledPoint> &scaled_points,
    const std::vector<size_t> &group_ids, std::int64_t scale) {
  EnsureNonEmptyPoints(scaled_points);
  ValidateOverflowFreeExpansion(scaled_points, scale);
  if (group_ids.size() != scaled_points.size()) {
    throw std::invalid_argument("group_ids length must match number of points");
  }

  std::vector<std::vector<size_t>> partition;
  std::unordered_map<size_t, size_t> group_index;
  partition.reserve(group_ids.size());

  for (size_t idx = 0; idx < group_ids.size(); ++idx) {
    const size_t group_id = group_ids[idx];
    auto [it, inserted] = group_index.emplace(group_id, partition.size());
    if (inserted) {
      partition.push_back({});
    }
    partition[it->second].push_back(idx);
  }
  auto result = SearchPartition(scaled_points, partition, scale);
  CoverStatus status = CoverStatus::kNonCoverable;
  switch (result.status) {
  case SearchStatus::kCoverable:
    status = CoverStatus::kCoverable;
    logging::Log(logging::Level::kInfo,
                 "Found disjoint cover for provided groups");
    break;
  case SearchStatus::kNonCoverable:
    logging::Log(logging::Level::kInfo,
                 "No disjoint cover found for provided groups");
    break;
  case SearchStatus::kUndecided:
    status = CoverStatus::kUndecided;
    logging::Log(logging::Level::kInfo, "Search undecided for provided groups");
    break;
  }

  return DisjointCoverReport{
      .status = status,
      .scale = scale,
      .points_scaled = scaled_points,
      .partitions_evaluated = std::nullopt,
      .sample_partition = result.status == SearchStatus::kCoverable
                              ? std::optional{partition}
                              : std::nullopt,
      .sample_centers = result.sample_centers,
  };
}

} // namespace disjoint_unit_disk
