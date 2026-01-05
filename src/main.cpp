#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>
#include <toml++/toml.h>

#include "disjoint_unit_disk/cli_io.hpp"
#include "disjoint_unit_disk/int_scaled_point.hpp"
#include "disjoint_unit_disk/logging.hpp"
#include "disjoint_unit_disk/search.hpp"
#include "disjoint_unit_disk/workflow.hpp"

namespace {

std::string JoinIndices(const std::vector<size_t> &indices) {
  std::string result = "(";
  for (size_t i = 0; i < indices.size(); ++i) {
    result += std::to_string(indices[i]);
    if (i + 1 < indices.size()) {
      result += ",";
    }
  }
  result += ")";
  return result;
}

void PrintReportSummary(const disjoint_unit_disk::DisjointCoverReport &report) {
  const double scale_d = static_cast<double>(report.scale);
  std::cout << "Points analyzed (" << report.points_scaled.size() << ")"
            << std::endl;
  for (size_t idx = 0; idx < report.points_scaled.size(); ++idx) {
    const auto &pt = report.points_scaled[idx];
    std::cout << "  " << idx << ": (" << static_cast<double>(pt.x) / scale_d
              << ", " << static_cast<double>(pt.y) / scale_d << ")"
              << std::endl;
  }

  if (!report.sample_partition.has_value() ||
      !report.sample_centers.has_value()) {
    switch (report.status) {
    case disjoint_unit_disk::CoverStatus::kNonCoverable:
      std::cout << "No disjoint disk cover found for the provided points."
                << std::endl;
      break;
    case disjoint_unit_disk::CoverStatus::kUndecided:
      std::cout << "Search was undecided; no sample cover available."
                << std::endl;
      break;
    case disjoint_unit_disk::CoverStatus::kCoverable:
    default:
      std::cout << "No sample disjoint disk cover is available." << std::endl;
      break;
    }
    return;
  }

  const auto &partition = *report.sample_partition;
  const auto &centers = report.sample_centers->centers;
  const double centers_scale_d =
      static_cast<double>(report.sample_centers->scale);
  std::cout << "Sample disjoint disk cover:" << std::endl;
  const auto precise_digits = std::numeric_limits<double>::max_digits10;
  std::cout << std::setprecision(precise_digits);
  for (size_t idx = 0; idx < partition.size(); ++idx) {
    const auto &group = partition[idx];
    const auto &center = centers[idx];
    std::cout << "  Disk " << idx << " covers " << JoinIndices(group)
              << " center=(" << static_cast<double>(center.x) / centers_scale_d
              << ", " << static_cast<double>(center.y) / centers_scale_d << ")"
              << std::endl;
  }
}

} // namespace

toml::table BuildReportToml(const disjoint_unit_disk::DisjointCoverReport &report,
                            const disjoint_unit_disk::ParsedPoints &parsed) {
  toml::table root;

  auto status_to_string = [](disjoint_unit_disk::CoverStatus s) -> std::string_view {
    switch (s) {
    case disjoint_unit_disk::CoverStatus::kCoverable:
      return "coverable";
    case disjoint_unit_disk::CoverStatus::kUndecided:
      return "undecided";
    case disjoint_unit_disk::CoverStatus::kNonCoverable:
    default:
      return "non_coverable";
    }
  };

  root.insert("status", status_to_string(report.status));
  if (report.partitions_evaluated.has_value()) {
    root.insert("partitions_evaluated",
                static_cast<std::int64_t>(*report.partitions_evaluated));
  }

  toml::table input_tbl;
  input_tbl.insert("scale", report.scale);
  toml::array points_arr;
  for (size_t i = 0; i < report.points_scaled.size(); ++i) {
    const auto &pt = report.points_scaled[i];
    toml::table pt_tbl;
    pt_tbl.insert("x", pt.x);
    pt_tbl.insert("y", pt.y);
    if (parsed.group_ids && i < parsed.group_ids->size()) {
      pt_tbl.insert("group",
                    static_cast<std::int64_t>((*parsed.group_ids)[i]));
    }
    points_arr.push_back(std::move(pt_tbl));
  }
  input_tbl.insert("points", std::move(points_arr));
  root.insert("input", std::move(input_tbl));

  if (report.sample_partition && report.sample_centers) {
    toml::table sample_tbl;
    toml::array partition_arr;
    for (const auto &group : *report.sample_partition) {
      toml::table part_tbl;
      toml::array indices;
      for (size_t idx : group) {
        indices.push_back(static_cast<std::int64_t>(idx));
      }
      part_tbl.insert("points", std::move(indices));
      partition_arr.push_back(std::move(part_tbl));
    }
    sample_tbl.insert("partition", std::move(partition_arr));

    toml::array centers_arr;
    for (const auto &center : report.sample_centers->centers) {
      toml::table c_tbl;
      c_tbl.insert("x", center.x);
      c_tbl.insert("y", center.y);
      centers_arr.push_back(std::move(c_tbl));
    }
    sample_tbl.insert("centers", std::move(centers_arr));
    sample_tbl.insert("scale", report.sample_centers->scale);

    root.insert("sample", std::move(sample_tbl));
  }

  return root;
}

int main(int argc, char **argv) {
  std::unique_ptr<CLI::App> app_ptr;
  try {
    app_ptr = std::make_unique<CLI::App>(
        "Analyze whether a planar point set can be covered by disjoint "
        "unit-radius disks.");
    auto &app = *app_ptr;

    std::string points_path;
    app.add_option("--points-path", points_path,
                   "TOML file containing [[points]] entries with x, y, and "
                   "optional group.")
        ->type_name("PATH")
        ->required();

    std::string report_path;
    app.add_option("--report-path", report_path,
                   "If set, write analysis results to a TOML file at this path.")
        ->type_name("PATH");

    bool verbose = false;
    app.add_flag("--verbose", verbose,
                 "Enable info-level logging during enumeration/search.");

    app.parse(argc, argv);

    disjoint_unit_disk::ParsedPoints parsed;
    try {
      parsed = disjoint_unit_disk::ParsePointsFile(points_path);
    } catch (const std::exception &ex) {
      std::cerr << "Error: " << ex.what() << std::endl;
      return 64; // EX_USAGE
    }

    disjoint_unit_disk::logging::SetMinimumLevel(
        verbose ? disjoint_unit_disk::logging::Level::kInfo
                : disjoint_unit_disk::logging::Level::kWarning);

    disjoint_unit_disk::DisjointCoverReport report;
    if (parsed.group_ids.has_value()) {
      report = disjoint_unit_disk::AnalyzeDisjointUnitDiskCoverWithGroups(
          parsed.points, *parsed.group_ids, parsed.scale);
    } else {
      report = disjoint_unit_disk::AnalyzeDisjointUnitDiskCover(parsed.points,
                                                           parsed.scale);
    }

    PrintReportSummary(report);
    if (!report_path.empty()) {
      try {
        toml::table report_toml = BuildReportToml(report, parsed);
        std::ofstream ofs(report_path, std::ios::trunc);
        if (!ofs) {
          std::cerr << "Error: failed to open report path: " << report_path
                    << std::endl;
          return 70; // EX_SOFTWARE
        }
        ofs << report_toml;
      } catch (const std::exception &ex) {
        std::cerr << "Error: failed to write report: " << ex.what()
                  << std::endl;
        return 70; // EX_SOFTWARE
      }
    }
    return report.status == disjoint_unit_disk::CoverStatus::kCoverable ? 0 : 1;
  } catch (const CLI::ParseError &e) {
    if (app_ptr) {
      return app_ptr->exit(e);
    }
    return 64; // EX_USAGE fallback if app construction failed
  } catch (const std::exception &ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
    return 70; // EX_SOFTWARE
  }
}
