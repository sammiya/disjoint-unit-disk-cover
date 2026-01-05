#include <gtest/gtest.h>

#include "disjoint_unit_disk/int_scaled_point.hpp"

namespace disjoint_unit_disk {
namespace {

// Minimal sanity: int distance helper
TEST(PointTest, DistanceSquaredIntPoint) {
  IntScaledPoint a{3, 4};
  IntScaledPoint b{-1, 1};
  EXPECT_EQ(IntDistanceSquared(a, b), 25);
}

} // namespace
} // namespace disjoint_unit_disk
