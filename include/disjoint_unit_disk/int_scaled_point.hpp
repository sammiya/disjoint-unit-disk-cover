#pragma once

#include <cstdint>

namespace disjoint_unit_disk {

struct IntScaledPoint {
  std::int64_t x;
  std::int64_t y;
};

inline std::int64_t IntDistanceSquared(const IntScaledPoint &a,
                                       const IntScaledPoint &b) {
  const std::int64_t dx = a.x - b.x;
  const std::int64_t dy = a.y - b.y;
  return dx * dx + dy * dy;
}

} // namespace disjoint_unit_disk
