#include "disjoint_unit_disk/cli_io.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include <toml++/toml.h>

namespace disjoint_unit_disk {
namespace {

ParsedPoints ParsePointsFileImpl(const std::string &file_path) {
  const std::filesystem::path path = file_path;
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("Points file not found: " + path.string());
  }

  toml::table data;
  try {
    data = toml::parse_file(path.string());
  } catch (const toml::parse_error &err) {
    throw std::runtime_error("Failed to parse TOML file '" + path.string() +
                             "': " + std::string(err.description()));
  }

  auto scale_opt = data["scale"].value<int64_t>();
  if (!scale_opt) {
    throw std::runtime_error(path.string() +
                             ": missing required integer field 'scale'");
  }
  if (*scale_opt < 1) {
    throw std::runtime_error(path.string() + ": scale must be an integer >= 1");
  }

  const toml::node_view points_node = data["points"];
  if (!points_node || !points_node.is_array()) {
    throw std::runtime_error(path.string() +
                             ": expected an array of [[points]]");
  }

  ParsedPoints parsed;
  parsed.scale = *scale_opt;
  bool saw_group = false;
  bool saw_no_group = false;

  const auto *arr = points_node.as_array();
  for (size_t idx = 0; idx < arr->size(); ++idx) {
    const auto *entry = (*arr)[idx].as_table();
    if (!entry) {
      throw std::runtime_error(path.string() + ": points[" +
                               std::to_string(idx) + "] is not a table");
    }

    auto x_opt = (*entry)["x"].value<int64_t>();
    auto y_opt = (*entry)["y"].value<int64_t>();
    if (!x_opt || !y_opt) {
      throw std::runtime_error(path.string() + ": points[" +
                               std::to_string(idx) + "] missing x or y");
    }
    parsed.points.push_back(IntScaledPoint{*x_opt, *y_opt});

    if (auto g_opt = (*entry)["group"].value<int64_t>()) {
      if (*g_opt < 0) {
        throw std::runtime_error(path.string() + ": points[" +
                                 std::to_string(idx) +
                                 "] group must be non-negative");
      }
      saw_group = true;
      if (!parsed.group_ids) {
        parsed.group_ids = std::vector<size_t>{};
        parsed.group_ids->reserve(arr->size());
      }
      parsed.group_ids->push_back(static_cast<size_t>(*g_opt));
    } else {
      saw_no_group = true;
    }
  }

  if (parsed.points.empty()) {
    throw std::runtime_error(path.string() +
                             ": file did not contain any points");
  }
  if (saw_group && saw_no_group) {
    throw std::runtime_error(path.string() +
                             ": group must be provided for all points or none");
  }
  if (saw_group && parsed.group_ids &&
      parsed.group_ids->size() != parsed.points.size()) {
    throw std::runtime_error(path.string() +
                             ": internal inconsistency parsing groups");
  }
  return parsed;
}

} // namespace

ParsedPoints ParsePointsFile(const std::string &points_file) {
  if (points_file.empty()) {
    throw std::runtime_error("No points file provided.");
  }
  return ParsePointsFileImpl(points_file);
}

} // namespace disjoint_unit_disk
