#include <gtest/gtest.h>

#include <vector>

#include "disjoint_unit_disk/search.hpp"

namespace disjoint_unit_disk {
namespace {

TEST(SearchTest, ReturnsDisksForSeparatedSingletonGroups) {
  std::vector<std::vector<IntScaledPoint>> groups{{IntScaledPoint{-1, 0}},
                                                  {IntScaledPoint{1, 0}}};
  auto result = Search(groups, 10);
  ASSERT_EQ(result.status, SearchStatus::kCoverable);
  ASSERT_TRUE(result.sample_centers.has_value());
  ASSERT_EQ(result.sample_centers->centers.size(), groups.size());
  const auto &centers = result.sample_centers->centers;
  const double scale_d = static_cast<double>(result.sample_centers->scale);
  const double dx = static_cast<double>(centers[0].x - centers[1].x) / scale_d;
  const double dy = static_cast<double>(centers[0].y - centers[1].y) / scale_d;
  const double distance_sq = dx * dx + dy * dy;
  EXPECT_GE(distance_sq, 4.0);
}

TEST(SearchTest, RectanglesPruningWhenGroupsTooClose) {
  std::vector<IntScaledPoint> group_a{
      {0, 0},
      {0, 19},
      {19, 0},
      {19, 19},
  };
  std::vector<IntScaledPoint> group_b{
      {8, 0},
      {8, 19},
      {27, 0},
      {27, 19},
  };
  EXPECT_EQ(Search({group_a, group_b}, 10).status, SearchStatus::kNonCoverable);
}

TEST(SearchTest, ReturnsNoneWhenGroupRectanglesDisjoint) {
  std::vector<IntScaledPoint> group{{0, 0}, {50, 0}};
  EXPECT_EQ(Search({group}, 10).status, SearchStatus::kNonCoverable);
}

TEST(SearchTest, RequiresAtLeastOneGroup) {
  EXPECT_THROW(Search({}, 1), std::invalid_argument);
}

TEST(SearchTest, RejectsEmptyGroup) {
  EXPECT_THROW(Search({{IntScaledPoint{0, 0}}, {}}, 1), std::invalid_argument);
}

TEST(SearchTest, NonIntegerCenterInstanceReturnsUndecided) {
  // Three points (8,0), (0,8), (11,11) are coverable by a radius-7 open disk
  // only when the center is non-lattice. The integer-center search reports
  // kUndecided; no sample centers are returned.
  std::vector<std::vector<IntScaledPoint>> groups{
      {{8, 0}, {0, 8}, {11, 11}},
  };
  const auto result = Search(groups, 7);
  EXPECT_EQ(result.status, SearchStatus::kUndecided);
  ASSERT_TRUE(result.sample_centers.has_value());
  EXPECT_TRUE(result.sample_centers->centers.empty());
}

} // namespace
} // namespace disjoint_unit_disk
