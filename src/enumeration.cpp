#include "disjoint_unit_disk/enumeration.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>

#include <CGAL/Exact_circular_kernel_2.h>

namespace disjoint_unit_disk {

namespace {

using CK = CGAL::Exact_circular_kernel_2;
using CircPoint = CK::Point_2;
using CircCircle = CGAL::Circle_2<CK>;
using ArcPoint = CGAL::Circular_arc_point_2<CK>;
using IntersectionResult =
    CGAL::CK2_Intersection_traits<CK, CircCircle, CircCircle>::type;

void BacktrackCoverPartitions(
    size_t n, const std::vector<std::vector<size_t>> &cover_sets,
    const std::vector<std::vector<size_t>> &subsets_by_element,
    std::vector<bool> &used, size_t covered,
    std::vector<std::vector<size_t>> &partition,
    std::vector<std::vector<std::vector<size_t>>> &results) {
  if (covered == n) {
    results.push_back(partition);
    return;
  }

  size_t next_index = n;
  for (size_t idx = 0; idx < n; ++idx) {
    if (!used[idx]) {
      next_index = idx;
      break;
    }
  }
  if (next_index == n) {
    return;
  }

  for (size_t subset_idx : subsets_by_element[next_index]) {
    const std::vector<size_t> &subset = cover_sets[subset_idx];
    bool conflict = false;
    for (size_t member : subset) {
      if (used[member]) {
        conflict = true;
        break;
      }
    }
    if (conflict) {
      continue;
    }

    for (size_t member : subset) {
      used[member] = true;
    }
    partition.push_back(subset);
    BacktrackCoverPartitions(n, cover_sets, subsets_by_element, used,
                             covered + subset.size(), partition, results);
    partition.pop_back();
    for (size_t member : subset) {
      used[member] = false;
    }
  }
}

} // namespace

std::vector<std::vector<size_t>>
EnumerateDiskCovers(const std::vector<IntScaledPoint> &scaled_points,
                    std::int64_t scale) {
  const size_t n = scaled_points.size();
  std::set<std::vector<size_t>> covers;
  const CK::FT scale_ft(scale);
  const CK::FT scale_sq(scale_ft * scale_ft);
  const std::int64_t scale_sq_int = scale * scale;

  // Lift input points into the circular kernel once; the original double values
  // are preserved as exact rationals.
  std::vector<CircPoint> cpoints;
  cpoints.reserve(n);
  for (const auto &p : scaled_points) {
    cpoints.emplace_back(CK::FT(p.x), CK::FT(p.y));
  }

  const auto squared_distance_int = [&](size_t a, size_t b) {
    const std::int64_t dx = scaled_points[a].x - scaled_points[b].x;
    const std::int64_t dy = scaled_points[a].y - scaled_points[b].y;
    return dx * dx + dy * dy;
  };

  auto emit_covers = [&](const std::vector<size_t> &interior,
                         const std::vector<size_t> &boundary) {
    const size_t boundary_count = boundary.size();
    const size_t max_bits = std::numeric_limits<size_t>::digits;
    if (boundary_count >= max_bits) {
      throw std::runtime_error(
          "Too many boundary points to enumerate subsets: boundary_count=" +
          std::to_string(boundary_count) +
          ", limit=" + std::to_string(max_bits - 1));
    }
    const size_t subset_total = static_cast<size_t>(1) << boundary_count;
    for (size_t mask = 0; mask < subset_total; ++mask) {
      std::vector<size_t> cover = interior;
      cover.reserve(interior.size() + boundary_count);
      for (size_t b = 0; b < boundary_count; ++b) {
        if (mask & (static_cast<size_t>(1) << b)) {
          cover.push_back(boundary[b]);
        }
      }
      if (!cover.empty()) {
        std::sort(cover.begin(), cover.end());
        covers.insert(std::move(cover));
      }
    }
  };

  // Centers placed at the original points (exact classification in CK).
  for (size_t i = 0; i < n; ++i) {
    std::vector<size_t> interior;
    std::vector<size_t> boundary;

    for (size_t j = 0; j < n; ++j) {
      const auto dist_sq = squared_distance_int(i, j);
      if (dist_sq < scale_sq_int) {
        interior.push_back(j);
      } else if (dist_sq == scale_sq_int) {
        boundary.push_back(j);
      }
    }

    emit_covers(interior, boundary);
  }

  // Centers at circle intersections of point pairs (all computations in the
  // circular kernel).
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      const auto dist_sq = squared_distance_int(i, j);
      if (dist_sq > 4 * scale_sq_int) {
        continue;
      }

      CircCircle circle_i(cpoints[i], scale_sq);
      CircCircle circle_j(cpoints[j], scale_sq);
      std::vector<IntersectionResult> inters;
      CGAL::intersection(circle_i, circle_j, std::back_inserter(inters));

      using PairAP = std::pair<ArcPoint, unsigned>;
      for (const auto &obj : inters) {
        if (const auto *ap = std::get_if<PairAP>(&obj)) {
          const auto &center_ap = ap->first;
          const auto cx = center_ap.x();
          const auto cy = center_ap.y();

          std::vector<size_t> interior;
          std::vector<size_t> boundary;
          interior.reserve(n);
          boundary.reserve(n);

          for (size_t k = 0; k < n; ++k) {
            if (k == i || k == j) {
              continue;
            }
            const auto dx = cx - cpoints[k].x();
            const auto dy = cy - cpoints[k].y();
            const auto cmp = CGAL::compare(dx * dx + dy * dy, scale_sq);
            if (cmp == CGAL::SMALLER) {
              interior.push_back(k);
            } else if (cmp == CGAL::EQUAL) {
              boundary.push_back(k);
            }
          }

          // The two generating points always lie on the boundary.
          boundary.push_back(i);
          boundary.push_back(j);

          emit_covers(interior, boundary);
        }
      }
    }
  }

  return std::vector<std::vector<size_t>>(covers.begin(), covers.end());
}

std::vector<std::vector<std::vector<size_t>>>
EnumerateCoverPartitions(const std::vector<std::vector<size_t>> &cover_sets,
                         size_t n) {
  if (n == 0) {
    throw std::invalid_argument("n must be positive");
  }

  std::vector<std::vector<size_t>> subsets_by_element(n);
  for (size_t idx = 0; idx < cover_sets.size(); ++idx) {
    for (size_t element : cover_sets[idx]) {
      if (element >= n) {
        throw std::out_of_range("cover set index out of range");
      }
      subsets_by_element[element].push_back(idx);
    }
  }

  std::vector<std::vector<std::vector<size_t>>> results;
  std::vector<std::vector<size_t>> partition;
  std::vector<bool> used(n, false);
  size_t covered = 0;

  BacktrackCoverPartitions(n, cover_sets, subsets_by_element, used, covered,
                           partition, results);
  return results;
}

} // namespace disjoint_unit_disk
