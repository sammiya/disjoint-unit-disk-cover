#include "disjoint_unit_disk/cli_io.hpp"

#include <filesystem>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

using disjoint_unit_disk::ParsedPoints;
using disjoint_unit_disk::ParsePointsFile;

namespace {

std::string WriteTempToml(const std::string &contents) {
  const auto tmpdir = std::filesystem::temp_directory_path();
  auto path = tmpdir / std::filesystem::path("ddc_test_points.toml");
  std::ofstream ofs(path);
  ofs << contents;
  ofs.close();
  return path.string();
}

} // namespace

TEST(CliIoTest, FilePointsWithoutGroups) {
  const std::string file_path = WriteTempToml(R"(
scale = 1

[[points]]
x = 0
y = 1

[[points]]
x = 2
y = 3
)");
  ParsedPoints parsed = ParsePointsFile(file_path);
  ASSERT_EQ(parsed.points.size(), 2u);
  EXPECT_FALSE(parsed.group_ids.has_value());
}

TEST(CliIoTest, FilePointsWithGroups) {
  const std::string file_path = WriteTempToml(R"(
scale = 2

[[points]]
x = 0
y = 1
group = 2

[[points]]
x = 2
y = 3
group = 4
)");
  ParsedPoints parsed = ParsePointsFile(file_path);
  ASSERT_TRUE(parsed.group_ids.has_value());
  const auto &groups = parsed.group_ids.value();
  EXPECT_EQ(parsed.scale, 2);
  ASSERT_EQ(groups.size(), 2u);
  EXPECT_EQ(groups[0], 2u);
  EXPECT_EQ(groups[1], 4u);
}

TEST(CliIoTest, FileMixedGroupsThrows) {
  const std::string file_path = WriteTempToml(R"(
scale = 3

[[points]]
x = 0
y = 1
group = 2

[[points]]
x = 2
y = 3
)");
  EXPECT_THROW(ParsePointsFile(file_path), std::runtime_error);
}

TEST(CliIoTest, MissingFileThrows) {
  EXPECT_THROW(ParsePointsFile("/nonexistent/path/to/points.toml"),
               std::runtime_error);
}

TEST(CliIoTest, MissingScaleThrows) {
  const std::string file_path = WriteTempToml(R"(
[[points]]
x = 0
y = 1
)");
  EXPECT_THROW(ParsePointsFile(file_path), std::runtime_error);
}

TEST(CliIoTest, InvalidScaleThrows) {
  const std::string file_path = WriteTempToml(R"(
scale = 0

[[points]]
x = 0
y = 1
)");
  EXPECT_THROW(ParsePointsFile(file_path), std::runtime_error);
}
