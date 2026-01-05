#include "disjoint_unit_disk/search.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>

namespace disjoint_unit_disk {
namespace {

struct IntRect {
  std::int64_t xmin;
  std::int64_t ymin;
  std::int64_t xmax;
  std::int64_t ymax;

  std::int64_t width() const { return xmax - xmin; }
  std::int64_t height() const { return ymax - ymin; }
};

std::optional<IntRect> IntersectRectangles(const std::vector<IntRect> &rects) {
  if (rects.empty()) {
    throw std::invalid_argument(
        "at least one rectangle is required for intersection");
  }
  IntRect current = rects.front();
  for (size_t i = 1; i < rects.size(); ++i) {
    const auto &r = rects[i];
    const std::int64_t ixmin = std::max(current.xmin, r.xmin);
    const std::int64_t iymin = std::max(current.ymin, r.ymin);
    const std::int64_t ixmax = std::min(current.xmax, r.xmax);
    const std::int64_t iymax = std::min(current.ymax, r.ymax);
    if (!(ixmin < ixmax) || !(iymin < iymax)) {
      return std::nullopt;
    }
    current = {ixmin, iymin, ixmax, iymax};
  }
  return current;
}

std::optional<std::array<IntRect, 2>> SplitLongestSide(const IntRect &rect) {
  const bool split_horizontal = rect.width() >= rect.height();
  if (split_horizontal) {
    const std::int64_t mid = std::midpoint(rect.xmin, rect.xmax);
    if (mid == rect.xmin || mid == rect.xmax) {
      return std::nullopt;
    }
    return std::array<IntRect, 2>{
        IntRect{rect.xmin, rect.ymin, mid, rect.ymax},
        IntRect{mid, rect.ymin, rect.xmax, rect.ymax},
    };
  }

  const std::int64_t mid = std::midpoint(rect.ymin, rect.ymax);
  if (mid == rect.ymin || mid == rect.ymax) {
    return std::nullopt;
  }
  return std::array<IntRect, 2>{
      IntRect{rect.xmin, rect.ymin, rect.xmax, mid},
      IntRect{rect.xmin, mid, rect.xmax, rect.ymax},
  };
}

IntScaledPoint RectCenter(const IntRect &rect) {
  const std::int64_t cx = std::midpoint(rect.xmin, rect.xmax);
  const std::int64_t cy = std::midpoint(rect.ymin, rect.ymax);
  return IntScaledPoint{cx, cy};
}

bool PointRectDistanceLess(const IntRect &rect, const IntScaledPoint &point,
                           std::int64_t threshold_sq) {
  const std::int64_t clamped_x = std::clamp(point.x, rect.xmin, rect.xmax);
  const std::int64_t clamped_y = std::clamp(point.y, rect.ymin, rect.ymax);
  const IntScaledPoint closest{clamped_x, clamped_y};
  return IntDistanceSquared(closest, point) < threshold_sq;
}

bool RectanglesFarthestDistanceLessThan(const IntRect &a, const IntRect &b,
                                        std::int64_t threshold_sq) {
  const std::array<IntScaledPoint, 4> corners_a{
      IntScaledPoint{a.xmin, a.ymin},
      IntScaledPoint{a.xmin, a.ymax},
      IntScaledPoint{a.xmax, a.ymin},
      IntScaledPoint{a.xmax, a.ymax},
  };
  const std::array<IntScaledPoint, 4> corners_b{
      IntScaledPoint{b.xmin, b.ymin},
      IntScaledPoint{b.xmin, b.ymax},
      IntScaledPoint{b.xmax, b.ymin},
      IntScaledPoint{b.xmax, b.ymax},
  };

  for (const auto &pa : corners_a) {
    for (const auto &pb : corners_b) {
      if (IntDistanceSquared(pa, pb) >= threshold_sq) {
        return false;
      }
    }
  }
  return true;
}

struct CacheState {
  std::vector<std::vector<bool>> distance_ok;
  std::vector<std::vector<std::uint64_t>> distance_entry_versions;
  std::vector<std::uint64_t> distance_row_versions;
  std::vector<std::vector<bool>> coverage_values;
  std::vector<std::vector<std::uint64_t>> coverage_entry_versions;
  std::vector<std::uint64_t> coverage_row_versions;

  void InvalidateRow(size_t row) {
    ++distance_row_versions[row];
    ++coverage_row_versions[row];

    // Handling overflow of version counters
    if (distance_row_versions[row] == 0) {
      distance_row_versions[row] = 1;
      for (auto &v : distance_entry_versions[row]) {
        v = 0;
      }
    }
    if (coverage_row_versions[row] == 0) {
      coverage_row_versions[row] = 1;
      for (auto &v : coverage_entry_versions[row]) {
        v = 0;
      }
    }
  }
};

struct SearchProgress {
};

SearchStatus
SearchRectangles(const std::vector<std::vector<IntScaledPoint>> &groups,
                 std::vector<IntRect> &rects,
                 std::vector<std::int64_t> &longest_sides,
                 std::vector<IntScaledPoint> &centers, CacheState &cache,
                 std::optional<size_t> last_split_index, std::int64_t scale,
                 std::vector<IntScaledPoint> &solution) {
  const size_t n = groups.size();
  const std::int64_t scale_sq = scale * scale;
  const std::int64_t four_scale_sq = 4 * scale_sq;

  bool covers_all = true;
  for (size_t i = 0; i < n && covers_all; ++i) {
    const auto &group = groups[i];
    const std::uint64_t row_version = cache.coverage_row_versions[i];
    for (size_t j = 0; j < group.size(); ++j) {
      if (cache.coverage_entry_versions[i][j] != row_version) {
        cache.coverage_values[i][j] =
            IntDistanceSquared(centers[i], group[j]) < scale_sq;
        cache.coverage_entry_versions[i][j] = row_version;
      }
      if (!cache.coverage_values[i][j]) {
        covers_all = false;
        break;
      }
    }
  }

  if (covers_all) {
    bool centers_ok = true;
    for (size_t i = 0; i < n && centers_ok; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        if (IntDistanceSquared(centers[i], centers[j]) < four_scale_sq) {
          centers_ok = false;
          break;
        }
      }
    }
    if (centers_ok) {
      solution = centers;
      return SearchStatus::kCoverable;
    }
  }

  for (size_t i = 0; i < n; ++i) {
    const auto &group = groups[i];
    const std::uint64_t row_version = cache.distance_row_versions[i];
    for (size_t j = 0; j < group.size(); ++j) {
      if (cache.distance_entry_versions[i][j] != row_version) {
        cache.distance_ok[i][j] =
            PointRectDistanceLess(rects[i], group[j], scale_sq);
        cache.distance_entry_versions[i][j] = row_version;
      }
      if (!cache.distance_ok[i][j]) {
        return SearchStatus::kNonCoverable;
      }
    }
  }

  const auto rectangles_too_close = [&](size_t first, size_t second) {
    return RectanglesFarthestDistanceLessThan(rects[first], rects[second],
                                              four_scale_sq);
  };

  if (last_split_index.has_value()) {
    const size_t i = *last_split_index;
    for (size_t j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      if (rectangles_too_close(i, j)) {
        return SearchStatus::kNonCoverable;
      }
    }
  } else {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        if (rectangles_too_close(i, j)) {
          return SearchStatus::kNonCoverable;
        }
      }
    }
  }

  const auto max_iter =
      std::max_element(longest_sides.begin(), longest_sides.end());
  const size_t idx =
      static_cast<size_t>(std::distance(longest_sides.begin(), max_iter));

  const IntRect original_rect = rects[idx];
  const std::int64_t original_longest = longest_sides[idx];
  const IntScaledPoint original_center = centers[idx];

  const auto subrects = SplitLongestSide(original_rect);
  if (!subrects.has_value()) {
    return SearchStatus::kUndecided;
  }

  SearchStatus best_child_status = SearchStatus::kNonCoverable;
  for (const IntRect &subrect : *subrects) {
    rects[idx] = subrect;
    longest_sides[idx] = std::max(subrect.width(), subrect.height());
    centers[idx] = RectCenter(subrect);
    cache.InvalidateRow(idx);

    const auto child_status =
        SearchRectangles(
            groups, rects, longest_sides, centers, cache, idx, scale, solution);
    if (child_status == SearchStatus::kCoverable) {
      return SearchStatus::kCoverable;
    }
    if (child_status == SearchStatus::kUndecided) {
      best_child_status = SearchStatus::kUndecided;
    }
  }

  rects[idx] = original_rect;
  longest_sides[idx] = original_longest;
  centers[idx] = original_center;
  cache.InvalidateRow(idx);

  return best_child_status;
}

} // namespace

SearchResult
Search(const std::vector<std::vector<IntScaledPoint>> &point_groups_scaled,
       std::int64_t scale) {
  if (point_groups_scaled.empty()) {
    throw std::invalid_argument("at least one group of points is required");
  }

  for (size_t idx = 0; idx < point_groups_scaled.size(); ++idx) {
    if (point_groups_scaled[idx].empty()) {
      throw std::invalid_argument("each group must have at least one point");
    }
  }

  std::vector<IntRect> feasible_rects;
  feasible_rects.reserve(point_groups_scaled.size());
  for (const auto &group : point_groups_scaled) {
    std::vector<IntRect> bounds;
    bounds.reserve(group.size());
    for (const auto &point : group) {
      bounds.push_back(IntRect{
          point.x - scale,
          point.y - scale,
          point.x + scale,
          point.y + scale,
      });
    }
    auto rect = IntersectRectangles(bounds);
    if (!rect.has_value()) {
      return SearchResult{SearchStatus::kNonCoverable, {}};
    }
    feasible_rects.push_back(*rect);
  }

  const size_t n = point_groups_scaled.size();
  std::vector<std::int64_t> longest_sides;
  longest_sides.reserve(n);
  for (const auto &rect : feasible_rects) {
    longest_sides.push_back(std::max(rect.width(), rect.height()));
  }

  std::vector<IntScaledPoint> centers;
  centers.reserve(n);
  for (const auto &rect : feasible_rects) {
    centers.push_back(RectCenter(rect));
  }

  CacheState cache;
  cache.distance_ok.resize(n);
  cache.distance_entry_versions.resize(n);
  cache.distance_row_versions.assign(n, 1);
  cache.coverage_values.resize(n);
  cache.coverage_entry_versions.resize(n);
  cache.coverage_row_versions.assign(n, 1);
  for (size_t i = 0; i < n; ++i) {
    const size_t group_size = point_groups_scaled[i].size();
    cache.distance_ok[i].assign(group_size, true);
    cache.distance_entry_versions[i].assign(group_size, 0);
    cache.coverage_values[i].assign(group_size, false);
    cache.coverage_entry_versions[i].assign(group_size, 0);
  }

  std::vector<IntScaledPoint> solution;
  const auto status = SearchRectangles(point_groups_scaled, feasible_rects,
                                       longest_sides, centers, cache,
                                       std::nullopt, scale, solution);
  switch (status) {
  case SearchStatus::kCoverable:
    return SearchResult{SearchStatus::kCoverable,
                        ScaledCenters{scale, std::move(solution)}};
  case SearchStatus::kUndecided:
    return SearchResult{SearchStatus::kUndecided,
                        ScaledCenters{scale, std::move(solution)}};
  case SearchStatus::kNonCoverable:
  default:
    return SearchResult{SearchStatus::kNonCoverable, std::nullopt};
  }
}

} // namespace disjoint_unit_disk
