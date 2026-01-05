#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

#include "disjoint_unit_disk/enumeration.hpp"

namespace disjoint_unit_disk {
namespace {

TEST(EnumerateUnitDiskCoversTest, SingletonsForWellSpacedPoints) {
  std::vector<IntScaledPoint> points{{0, 0}, {3, 0}, {-3, 0}};
  const auto result = EnumerateDiskCovers(points, 1);
  const std::vector<std::vector<size_t>> expected{{0}, {1}, {2}};
  EXPECT_EQ(result, expected);
}

TEST(EnumerateUnitDiskCoversTest, BoundaryPointsGenerateAllSubsets) {
  // Center at (0,0) with four boundary points exactly at distance 1.
  std::vector<IntScaledPoint> points{
      {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1},
  };

  const auto covers = EnumerateDiskCovers(points, 1);

  size_t count_center = 0;
  for (const auto &cover : covers) {
    if (!cover.empty() && cover[0] == 0) {
      ++count_center;
    }
  }
  EXPECT_EQ(count_center, 16u);
}

TEST(EnumerateUnitDiskCoversTest, EquilateralPointsWithinUnitDisk) {
  const double sqrt3 = std::sqrt(3.0);
  const std::int64_t scale_header = 100000000;
  const auto round_to_int = [](double v) {
    return static_cast<std::int64_t>(std::llround(v));
  };
  std::vector<IntScaledPoint> points{
      {round_to_int(0.99 * 1.0 * scale_header), 0},
      {round_to_int(0.99 * -0.5 * scale_header),
       round_to_int(0.99 * (sqrt3 / 2.0) * scale_header)},
      {round_to_int(0.99 * -0.5 * scale_header),
       round_to_int(0.99 * -(sqrt3 / 2.0) * scale_header)},
  };
  const auto result = EnumerateDiskCovers(points, scale_header);
  std::set<std::vector<size_t>> canonical(result.begin(), result.end());
  std::set<std::vector<size_t>> expected{{0}, {0, 1}, {0, 1, 2}, {0, 2},
                                         {1}, {1, 2}, {2}};
  EXPECT_EQ(canonical, expected);
}

TEST(EnumerateUnitDiskCoversTest, EquilateralPointsOutsideUnitDisk) {
  const double sqrt3 = std::sqrt(3.0);
  const std::int64_t scale_header = 100000000;
  const auto round_to_int = [](double v) {
    return static_cast<std::int64_t>(std::llround(v));
  };
  std::vector<IntScaledPoint> points{
      {round_to_int(1.01 * 1.0 * scale_header), 0},
      {round_to_int(1.01 * -0.5 * scale_header),
       round_to_int(1.01 * (sqrt3 / 2.0) * scale_header)},
      {round_to_int(1.01 * -0.5 * scale_header),
       round_to_int(1.01 * -(sqrt3 / 2.0) * scale_header)},
  };
  const auto result = EnumerateDiskCovers(points, scale_header);
  std::set<std::vector<size_t>> canonical(result.begin(), result.end());
  std::set<std::vector<size_t>> expected{
      {0}, {1}, {2}, {0, 1}, {0, 2}, {1, 2},
  };
  EXPECT_EQ(canonical, expected);
}

TEST(EnumerateUnitDiskCoversTest,
     TwoPointsWithinTwoUnitsProduceSingletonsAndPair) {
  std::vector<IntScaledPoint> points{
      {0, 0}, {3, 0}}; // distance 3 with scale=2 -> normalized 1.5
  const auto result = EnumerateDiskCovers(points, 2);

  std::set<std::vector<size_t>> canonical(result.begin(), result.end());
  std::set<std::vector<size_t>> expected{{0}, {1}, {0, 1}};
  EXPECT_EQ(canonical, expected);
}

TEST(EnumerateUnitDiskCoversTest, IsolatedPointAlwaysProducesSingleton) {
  // Point 2 is isolated: farther than 2 units from the others.
  std::vector<IntScaledPoint> points{
      {0, 0}, {3, 0}, {8, 0}}; // scale=2 -> normalized 0,1.5,4
  const auto result = EnumerateDiskCovers(points, 2);

  std::set<std::vector<size_t>> canonical(result.begin(), result.end());
  EXPECT_TRUE(canonical.count({2}) > 0);
}

TEST(EnumerateCoverPartitionsTest, PairwiseCandidatesFormPerfectMatchings) {
  std::vector<std::vector<size_t>> candidates{{0, 1}, {0, 2}, {0, 3},
                                              {1, 2}, {1, 3}, {2, 3}};
  const auto partitions = EnumerateCoverPartitions(candidates, 4);

  std::set<std::vector<std::vector<size_t>>> canonical;
  for (const auto &partition : partitions) {
    auto sorted_partition = partition;
    std::sort(sorted_partition.begin(), sorted_partition.end());
    canonical.insert(sorted_partition);
  }

  std::set<std::vector<std::vector<size_t>>> expected{
      {{0, 1}, {2, 3}},
      {{0, 2}, {1, 3}},
      {{0, 3}, {1, 2}},
  };

  EXPECT_EQ(canonical, expected);
}

TEST(EnumerateCoverPartitionsTest, NoPartitionWhenUncoveredIndex) {
  std::vector<std::vector<size_t>> candidates{{0}, {2}};
  const auto partitions = EnumerateCoverPartitions(candidates, 3);
  EXPECT_TRUE(partitions.empty());
}

} // namespace
} // namespace disjoint_unit_disk
