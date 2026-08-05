#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <regex>
#include <yaml-cpp/yaml.h>

#include "types.hpp"

/// @namespace ConfigUtils
/// @brief Utility functions for parsing configuration files
namespace ConfigUtils {

    /// @brief Parses measurement metadata from the provided configuration file path
    /// @param config_file_path The path to the YAML configuration file containing
    /// measurement metadata
    /// @param summary_root_file_out Optional pointer to a string where the parsed
    /// summary root file path will be stored
    /// @return A vector of parsed measurement metadata entries from config file
    inline std::vector<MeasurementMetadata> parseMeasurementMetadata(
            const std::string& config_file_path,
            std::string* summary_root_file_out = nullptr)
    {
        std::vector<MeasurementMetadata> metadata_entries;
        YAML::Node config = YAML::LoadFile(config_file_path);

        // Helper lambdas to read different types of values from YAML nodes
        auto read_bool = [](const YAML::Node& node, const char* key) -> bool {
            if (!node[key]) return false;
            if (node[key].IsScalar()) {
                std::string value = node[key].Scalar();
                return value == "true" || value == "1" || value == "ON";
            }
            return false;
        };
        auto read_string = [](const YAML::Node& node, const char* key) -> std::string {
            return (node[key] && node[key].IsScalar()) ? node[key].Scalar() : "";
        };
        auto read_double = [](const YAML::Node& node, const char* key) -> double {
            if (!node[key] || !node[key].IsScalar()) return 0.0;
            try { return std::stod(node[key].Scalar()); } catch (...) { return 0.0; }
        };
        auto extract_int = [](const std::string& str, int default_val = 0) -> int {
            std::smatch match;
            if (std::regex_search(str, match, std::regex("\\d+"))) {
                return std::stoi(match.str());
            }
            return default_val;
        };
        auto extract_float = [](const std::string& str, double default_val = 0.0) -> double {
            std::smatch match;
            if (std::regex_search(str, match, std::regex("\\d+([._]\\d+)?"))) {
                std::string num_str = match.str();
                std::replace(num_str.begin(), num_str.end(), '_', '.');
                try {
                    return std::stod(num_str);
                } catch (...) {
                    return default_val;
                }
            }
            return default_val;
        };

        // Read summary root file path from config, or construct default path if not provided
        std::string summary_root_file = read_string(config, "summary root file");
        if (summary_root_file.empty()) {
            std::filesystem::path config_stem = std::filesystem::path(config_file_path).stem();
            summary_root_file = (std::filesystem::path("data/output") / (config_stem.string() + "_summary.root")).string();
        }
        if (summary_root_file_out) *summary_root_file_out = summary_root_file;

        // Read global fields as STRINGS to check for the "SCAN" keyword
        std::string measurement_type  = read_string(config, "measurement_type");
        bool global_source            = read_bool(config, "source");

        std::string global_mixture    = read_string(config, "mixture");
        std::string global_lv         = read_string(config, "lv_setting");
        std::string global_filter_str = read_string(config, "filter");
        std::string global_layer      = read_string(config, "layer");

        YAML::Node data = config["data"];
        if (!data) return metadata_entries;

        // Helper lambda to build a MeasurementMetadata entry
        auto build_entry = [&](const YAML::Node& node, const std::string& group_key) {
            MeasurementMetadata entry;
            entry.name = read_string(node, "name");
            entry.measurement_type = measurement_type;
            entry.root_file = read_string(node, "root file");
            entry.group_name = group_key;

            entry.source = global_source;

            // RESOLVE MIXTURE
            if (global_mixture == "SCAN" || global_mixture == "scan") {
                entry.mixture = group_key; // e.g., "std_mixture"
            } else {
                entry.mixture = global_mixture;
            }

            // RESOLVE LV SETTING
            if (global_lv == "SCAN" || global_lv == "scan") {
                entry.lv_setting = extract_int(group_key); // Extracts 3 from "lv_3"
            } else {
                entry.lv_setting = extract_int(global_lv);
            }

            // RESOLVE FILTER
            if (global_filter_str == "SCAN" || global_filter_str == "scan") {
                if (group_key == "OFF" || group_key == "off" || group_key == "NONE" || group_key.find("filter_OFF") != std::string::npos) {
                    entry.filter = 0.0;
                } else {
                    entry.filter = extract_float(group_key); // Converts "filter_6_6" to 6.6
                }
            } else {
                if (global_filter_str == "OFF" || global_filter_str == "off" || global_filter_str == "NONE" || global_filter_str.empty()) {
                    entry.filter = 0.0;
                } else {
                    try { entry.filter = std::stod(global_filter_str); } catch (...) { entry.filter = -1.0; }
                }
            }

            // RESOLVE LAYER
            if (global_layer == "SCAN" || global_layer == "scan") {
                entry.scanned_layer = extract_int(group_key); // Extracts 0 from "layer0"
            } else {
                entry.scanned_layer = extract_int(global_layer);
            }

            // Point-specific values
            entry.scanned_hv = read_double(node, "hv");
            entry.other_hv = read_double(node, "other_hv");
            if (entry.scanned_hv == 0.0) {
                entry.scanned_hv = read_double(node, "scanned_hv");
                entry.other_hv = entry.scanned_hv;
            }

            if (!entry.name.empty()) metadata_entries.push_back(std::move(entry));
        };

        if (data.IsSequence()) {
            // Sequence format: no groupings, just a flat list of scans
            for (const auto& node : data) {
                build_entry(node, "default");
            }
        } else if (data.IsMap()) { 
            // Map format: data is grouped by keys
            for (const auto& item : data) {
                std::string group_key = item.first.Scalar();
                const YAML::Node& runs = item.second;
                if (runs.IsSequence()) {
                    for (const auto& node : runs) {
                        build_entry(node, group_key);
                    }
                }
            }
        }
        
        return metadata_entries;
    }
} // namespace ConfigUtils