#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <utility>
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

TEST(SearchTest, CenterDischargeHandlesSeveralGroups) {
  std::vector<std::vector<IntScaledPoint>> groups{
      {{-20, 0}}, {{0, 0}}, {{20, 0}},
  };
  for (bool center_discharge : {false, true}) {
    SCOPED_TRACE(center_discharge);
    const auto result = Search(groups, 10, center_discharge);
    ASSERT_EQ(result.status, SearchStatus::kCoverable);
    CheckWitness(groups, 10, result);
  }
}

TEST(SearchTest, CenterDischargePreservesUndecidedForNonIntegerCenter) {
  // Fixing a distant group's center must not change the remaining answer.
  std::vector<std::vector<IntScaledPoint>> groups{
      {{100, 100}}, {{8, 0}, {0, 8}, {11, 11}},
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

TEST(SearchTest, CenterDischargeDoesNotRejectSmallIntegerCovers) {
  // Exhaust all integer centers independently on reproducible small inputs.
  // A continuous cover with no integer witness may still be undecided.
  std::mt19937 random(20261003);
  constexpr std::int64_t scale = 3;
  for (int trial = 0; trial < 200; ++trial) {
    SCOPED_TRACE(trial);
    std::vector<std::vector<IntScaledPoint>> groups(3);
    std::vector<std::vector<IntScaledPoint>> possible_centers;
    for (auto &group : groups) {
      const auto x = static_cast<std::int64_t>(random() % 17) - 8;
      const auto y = static_cast<std::int64_t>(random() % 17) - 8;
      group = {{x, y}, {x + static_cast<std::int64_t>(random() % 5) - 2,
                       y + static_cast<std::int64_t>(random() % 5) - 2}};
      std::vector<IntScaledPoint> centers;
      for (auto cx = x - scale; cx <= x + scale; ++cx) {
        for (auto cy = y - scale; cy <= y + scale; ++cy) {
          const IntScaledPoint center{cx, cy};
          if (std::all_of(group.begin(), group.end(), [&](const auto &p) {
                return IntDistanceSquared(center, p) < scale * scale;
              })) {
            centers.push_back(center);
          }
        }
      }
      possible_centers.push_back(std::move(centers));
    }
    bool has_integer_cover = false;
    for (const auto &a : possible_centers[0]) {
      for (const auto &b : possible_centers[1]) {
        for (const auto &c : possible_centers[2]) {
          if (IntDistanceSquared(a, b) >= 4 * scale * scale &&
              IntDistanceSquared(a, c) >= 4 * scale * scale &&
              IntDistanceSquared(b, c) >= 4 * scale * scale) {
            has_integer_cover = true;
          }
        }
      }
    }
    for (bool center_discharge : {false, true}) {
      SCOPED_TRACE(center_discharge);
      const auto result = Search(groups, scale, center_discharge);
      if (has_integer_cover) {
        EXPECT_NE(result.status, SearchStatus::kNonCoverable);
      }
      if (result.status == SearchStatus::kCoverable) {
        EXPECT_TRUE(has_integer_cover);
        CheckWitness(groups, scale, result);
      }
    }
  }
}

} // namespace
} // namespace disjoint_unit_disk
