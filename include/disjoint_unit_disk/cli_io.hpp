#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "disjoint_unit_disk/int_scaled_point.hpp"

namespace disjoint_unit_disk {

struct ParsedPoints {
  std::int64_t scale;
  std::vector<IntScaledPoint> points;
  std::optional<std::vector<size_t>> group_ids;
};

ParsedPoints ParsePointsFile(const std::string &points_file);

} // namespace disjoint_unit_disk
