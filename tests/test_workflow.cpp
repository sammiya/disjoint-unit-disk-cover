#include <gtest/gtest.h>

#include <limits>
#include <vector>

#include "disjoint_unit_disk/int_scaled_point.hpp"
#include "disjoint_unit_disk/workflow.hpp"

namespace disjoint_unit_disk {
namespace {

constexpr std::int64_t kScale =
    10; // scales coordinates to integers (e.g., 0.6 -> 6)

std::vector<IntScaledPoint> BuildGridPointsScaled() {
  const int coords_scaled[] = {-6, 0, 6}; // coords * kScale
  std::vector<IntScaledPoint> points;
  for (int x : coords_scaled) {
    for (int y : coords_scaled) {
      points.push_back(IntScaledPoint{static_cast<std::int64_t>(x),
                                      static_cast<std::int64_t>(y)});
    }
  }
  points.push_back(IntScaledPoint{12, 6}); // (1.2, 0.6) scaled by 10
  return points;
}

TEST(WorkflowTest, AnalyzeDisjointCoverReturnsSample) {
  auto scaled_points = BuildGridPointsScaled();
  auto report = AnalyzeDisjointUnitDiskCover(scaled_points, kScale);
  ASSERT_EQ(report.status, CoverStatus::kCoverable);
  EXPECT_TRUE(report.sample_partition.has_value());
  EXPECT_TRUE(report.sample_centers.has_value());

  const auto &partition = *report.sample_partition;
  const auto &centers_scaled = report.sample_centers->centers;
  ASSERT_EQ(partition.size(), centers_scaled.size());

  const std::int64_t scale_sq =
      report.sample_centers->scale * report.sample_centers->scale;

  // Each center covers its assigned points within scale.
  for (size_t i = 0; i < partition.size(); ++i) {
    for (size_t pt_idx : partition[i]) {
      ASSERT_LT(pt_idx, scaled_points.size());
      EXPECT_LT(IntDistanceSquared(centers_scaled[i], scaled_points[pt_idx]),
                scale_sq);
    }
  }

  // Centers are pairwise at least 2*scale apart.
  const std::int64_t four_scale_sq = 4 * scale_sq;
  for (size_t i = 0; i + 1 < centers_scaled.size(); ++i) {
    for (size_t j = i + 1; j < centers_scaled.size(); ++j) {
      EXPECT_GE(IntDistanceSquared(centers_scaled[i], centers_scaled[j]),
                four_scale_sq);
    }
  }
}

TEST(WorkflowTest, AnalyzeDisjointCoverRequiresPoints) {
  EXPECT_THROW(AnalyzeDisjointUnitDiskCover({}, kScale), std::invalid_argument);
}

TEST(WorkflowTest, AnalyzeWithGroupsRunsSearch) {
  std::vector<IntScaledPoint> points{
      IntScaledPoint{-25, 0},
      IntScaledPoint{25, 0},
  };
  std::vector<size_t> groups{0, 1};
  auto report = AnalyzeDisjointUnitDiskCoverWithGroups(points, groups, kScale);
  ASSERT_EQ(report.status, CoverStatus::kCoverable);
  EXPECT_TRUE(report.sample_partition.has_value());
  EXPECT_TRUE(report.sample_centers.has_value());
}

TEST(WorkflowTest, AnalyzeWithGroupsRequiresMatchingLengths) {
  std::vector<IntScaledPoint> points{
      {0, 0},
      {10, 0},
  };
  EXPECT_THROW(AnalyzeDisjointUnitDiskCoverWithGroups(
                   points, std::vector<size_t>{0}, kScale),
               std::invalid_argument);
}

TEST(WorkflowTest, AnalyzeWithGroupsReturnsNoneWhenGroupsTooClose) {
  std::vector<IntScaledPoint> points{
      {0, 0}, {0, 19}, {19, 0}, {19, 19}, {8, 0}, {8, 19}, {27, 0}, {27, 19},
  };
  std::vector<size_t> group_ids{0, 0, 0, 0, 1, 1, 1, 1};

  auto report =
      AnalyzeDisjointUnitDiskCoverWithGroups(points, group_ids, kScale);
  EXPECT_EQ(report.status, CoverStatus::kNonCoverable);
  EXPECT_FALSE(report.sample_partition.has_value());
  EXPECT_FALSE(report.sample_centers.has_value());
}

TEST(WorkflowTest, AnalyzeWithGroupsReturnsUndecidedWhenOnlyNonIntegerCenter) {
  // Same instance as the search test: coverable by radius-7 with non-lattice
  // center only; integer-center search reports kUndecided with empty centers.
  std::vector<IntScaledPoint> points{{8, 0}, {0, 8}, {11, 11}};
  std::vector<size_t> group_ids{0, 0, 0};
  auto report =
      AnalyzeDisjointUnitDiskCoverWithGroups(points, group_ids, 7);
  EXPECT_EQ(report.status, CoverStatus::kUndecided);
  EXPECT_FALSE(report.sample_partition.has_value());
  ASSERT_TRUE(report.sample_centers.has_value());
  EXPECT_TRUE(report.sample_centers->centers.empty());
}

TEST(WorkflowTest, FailsWhenRectangleExpansionOverflows) {
  std::vector<IntScaledPoint> points{
      {std::numeric_limits<std::int64_t>::max(), 0},
  };
  EXPECT_THROW(AnalyzeDisjointUnitDiskCover(points, 1), std::runtime_error);
}

TEST(WorkflowTest, FailsWhenCornerDistancesOverflow) {
  const std::int64_t max = std::numeric_limits<std::int64_t>::max() - 1;
  const std::int64_t min = std::numeric_limits<std::int64_t>::min() + 1;
  std::vector<IntScaledPoint> points{
      {max, 0},
      {min, 0},
  };
  EXPECT_THROW(AnalyzeDisjointUnitDiskCover(points, 1), std::runtime_error);
}

} // namespace
} // namespace disjoint_unit_disk
