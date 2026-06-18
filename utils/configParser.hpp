#pragma once

#include <filesystem>
#include <string>
#include <vector>

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
        // Vector to hold parsed measurement metadata entries
        std::vector<MeasurementMetadata> metadata_entries;

        // Load the YAML configuration file
        YAML::Node config = YAML::LoadFile(config_file_path);

        // Helper lambdas to read different types of values from YAML nodes
        auto read_bool = [](const YAML::Node& node, const char* key) -> bool {
            if (!node[key]) {
                return false;
            }
            if (node[key].IsScalar()) {
                std::string value = node[key].Scalar();
                return value == "true" || value == "1" || value == "ON";
            }
            return false;
        };
        auto read_string = [](const YAML::Node& node, const char* key) -> std::string {
            if (!node[key]) {
                return "";
            }
            return node[key].IsScalar() ? node[key].Scalar() : "";
        };
        auto read_double = [](const YAML::Node& node, const char* key) -> double {
            if (!node[key] || !node[key].IsScalar()) {
                return 0.0;
            }
            try {
                return std::stod(node[key].Scalar());
            } catch (const std::exception&) {
                return 0.0;
            }
        };
        auto read_int = [](const YAML::Node& node, const char* key) -> int {
            if (!node[key] || !node[key].IsScalar()) {
                return 0;
            }
            try {
                return std::stoi(node[key].Scalar());
            } catch (const std::exception&) {
                return 0;
            }
        };
        auto read_filter = [](const YAML::Node& node, const char* key) -> double {
            // If the node doesn't exist or isn't a plain value, return default (OFF)
            if (!node[key] || !node[key].IsScalar()) {
                return 0.0; // Default to 0.0 for OFF
            }

            // Extract the raw string from YAML
            std::string value_str = node[key].Scalar();

            // Explicitly intercept the "OFF" token (case-sensitive check)
            if (value_str == "OFF" || value_str == "off" || value_str == "NONE") {
                return 0.0;
            }

            // If it isn't "OFF", try parsing it as a floating-point number
            try {
                return std::stod(value_str);
            } catch (const std::exception&) {
                // If parsing fails (e.g., a typo like "OOFF"), fall back to -1.0
                return -1.0; 
            }
        };

        // Read summary root file path from config, or construct default path if not provided
        std::string summary_root_file = read_string(config, "summary root file");
        if (summary_root_file.empty()) {
            std::filesystem::path config_stem = std::filesystem::path(config_file_path).stem();
            std::filesystem::path summary_root_path = std::filesystem::path("data/output") /
                                                    (config_stem.string() + "_summary.root");
            summary_root_file = summary_root_path.string();
        }
        if (summary_root_file_out) {
            *summary_root_file_out = summary_root_file;
        }

        // Read common metadata fields for all measurements in this config
        std::string measurement_type = read_string(config, "measurement_type");
        bool source = read_bool(config, "source");
        double filter = read_filter(config, "filter");
        std::string mixture = read_string(config, "mixture");
        int lv_setting = read_int(config, "lv_setting");
        YAML::Node data = config["data"];
        if (!data) return metadata_entries;

        // Helper lambda to build a MeasurementMetadata entry from a YAML node and layer number
        auto build_entry = [&](const YAML::Node& node, const int layer) {
            MeasurementMetadata entry;
            entry.name = read_string(node, "name");
            entry.measurement_type = measurement_type;
            entry.root_file = read_string(node, "root file");

            entry.mixture = mixture;
            entry.source = source;
            entry.filter = filter;
            entry.lv_setting = lv_setting;

            entry.scanned_layer = layer;
            entry.scanned_hv = read_double(node, "hv");
            entry.other_hv = read_double(node, "other_hv");

            if (entry.scanned_hv == 0.0) {
                entry.scanned_hv = read_double(node, "scanned_hv");
                entry.other_hv = entry.scanned_hv;
            }

            // Add the entry if it has a valid name (non-empty)
            if (!entry.name.empty()) metadata_entries.push_back(std::move(entry));
        };

        // Helper lambda to parse layer keys like "layer1", "layer2", etc. and extract the layer number
        auto parse_layer_key = [](const YAML::Node& key_node) -> int {
            if (key_node.IsScalar()) {
                std::string key = key_node.Scalar();
                if (key.rfind("layer", 0) == 0) {
                    key = key.substr(5);
                }
                try {
                    return std::stoi(key);
                } catch (const std::exception&) {
                    return 0;
                }
            }
            return 0;
        };

        if (data.IsSequence()) {   // Sequence format: data is a list of entries (noise scans)
            for (const auto& node : data) {
                build_entry(node, 1);
            }
        } else if (data.IsMap()) { // Map format: data is a map of layer keys to lists of entries (efficiency scans)
            for (const auto& item : data) {
                int layer = parse_layer_key(item.first);
                const YAML::Node& runs = item.second;
                if (runs.IsSequence()) {
                    for (const auto& node : runs) {
                        build_entry(node, layer);
                    }
                }
            }
        }
        return metadata_entries;
    }
} // namespace ConfigUtils