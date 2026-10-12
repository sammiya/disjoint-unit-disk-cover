#include <gtest/gtest.h>

#include <vector>

#include "disjoint_unit_disk/search.hpp"

namespace disjoint_unit_disk {
namespace {

void CheckWitness(const std::vector<std::vector<IntScaledPoint>> &groups,
                  std::int64_t scale, const SearchResult &result) {
  ASSERT_TRUE(result.sample_centers.has_value());
  const auto &centers = result.sample_centers->centers;
  ASSERT_EQ(centers.size(), groups.size());
  EXPECT_EQ(result.sample_centers->scale, scale);
  for (size_t i = 0; i < groups.size(); ++i) {
    for (const auto &point : groups[i]) {
      EXPECT_LT(IntDistanceSquared(centers[i], point), scale * scale);
    }
    for (size_t j = 0; j < i; ++j) {
      EXPECT_GE(IntDistanceSquared(centers[i], centers[j]), 4 * scale * scale);
    }
  }
}

TEST(SearchTest, ReturnsDisksForSeparatedSingletonGroups) {
  std::vector<std::vector<IntScaledPoint>> groups{{IntScaledPoint{-1, 0}},
                                                  {IntScaledPoint{1, 0}}};
  for (bool center_discharge : {false, true}) {
    SCOPED_TRACE(center_discharge);
    const auto result = Search(groups, 10, center_discharge);
    ASSERT_EQ(result.status, SearchStatus::kCoverable);
    CheckWitness(groups, 10, result);
  }
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
  for (bool center_discharge : {false, true}) {
    SCOPED_TRACE(center_discharge);
    EXPECT_EQ(Search({group_a, group_b}, 10, center_discharge).status,
              SearchStatus::kNonCoverable);
  }
}

TEST(SearchTest, ReturnsNoneWhenGroupRectanglesDisjoint) {
  std::vector<IntScaledPoint> group{{0, 0}, {50, 0}};
  for (bool center_discharge : {false, true}) {
    SCOPED_TRACE(center_discharge);
    EXPECT_EQ(Search({group}, 10, center_discharge).status,
              SearchStatus::kNonCoverable);
  }
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
  for (bool center_discharge : {false, true}) {
    SCOPED_TRACE(center_discharge);
    const auto result = Search(groups, 7, center_discharge);
    EXPECT_EQ(result.status, SearchStatus::kUndecided);
    ASSERT_TRUE(result.sample_centers.has_value());
    EXPECT_TRUE(result.sample_centers->centers.empty());
  }
}

TEST(SearchTest, CenterDischargeDoesNotCoverBoundaryPoints) {
  // Integer centers can cover the triangle only with boundary points.
  std::vector<std::vector<IntScaledPoint>> groups{
      {{100, 100}}, {{1, 0}, {0, 1}, {1, 1}},
  };
  for (bool center_discharge : {false, true}) {
    SCOPED_TRACE(center_discharge);
    const auto result = Search(groups, 1, center_discharge);
    EXPECT_EQ(result.status, SearchStatus::kUndecided);
    ASSERT_TRUE(result.sample_centers.has_value());
    EXPECT_TRUE(result.sample_centers->centers.empty());
  }
}

TEST(SearchTest, CenterDischargeHandlesInteractingGroups) {
  struct TestCase {
    const char *name;
    std::vector<std::vector<IntScaledPoint>> groups;
    SearchStatus expected_status;
  };
  const TestCase cases[]{
      // Integer centers have |y| <= 1, so disjoint disks require x-gaps of
      // at least 6, impossible for three centers within [-2,8].
      // Noninteger centers (-11/4,0), (3,7/4), (35/4,0) work.
      {"tight chain",
       {{{0, -1}, {0, 1}}, {{3, -1}, {3, 1}}, {{6, -1}, {6, 1}}},
       SearchStatus::kUndecided},
      // Centers (-2,-1), (2,4), (6,-1) give an integer cover.
      {"staggered groups",
       {{{0, -1}, {0, 1}}, {{2, 2}, {2, 4}}, {{4, -1}, {4, 1}}},
       SearchStatus::kCoverable},
  };
  constexpr std::int64_t scale = 3;
  for (const auto &test_case : cases) {
    SCOPED_TRACE(test_case.name);
    const auto &groups = test_case.groups;

    for (bool center_discharge : {false, true}) {
      SCOPED_TRACE(center_discharge);
      const auto result = Search(groups, scale, center_discharge);
      EXPECT_EQ(result.status, test_case.expected_status);
      if (result.status == SearchStatus::kCoverable) {
        CheckWitness(groups, scale, result);
      }
    }
  }
}

} // namespace
} // namespace disjoint_unit_disk
