#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "disjoint_unit_disk/int_scaled_point.hpp"

namespace disjoint_unit_disk {

std::vector<std::vector<size_t>>
EnumerateDiskCovers(const std::vector<IntScaledPoint> &scaled_points,
                    std::int64_t scale);

std::vector<std::vector<std::vector<size_t>>>
EnumerateCoverPartitions(const std::vector<std::vector<size_t>> &cover_sets,
                         size_t n);

} // namespace disjoint_unit_disk
