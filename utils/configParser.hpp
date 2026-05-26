#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "types.hpp"

/// @namespace ConfigUtils
/// @brief Utility functions for parsing configuration files
namespace ConfigUtils {
inline std::vector<MeasurementEntry> parseMeasurementEntries(const std::string& config_file_path,
                                                             std::string* summary_root_file_out = nullptr) {
    std::vector<MeasurementEntry> entries;

    YAML::Node config = YAML::LoadFile(config_file_path);

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

    std::string measurement_type = read_string(config, "measurement_type");
    std::string source = read_string(config, "source");
    double filter = read_double(config, "filter");
    std::string mixture = read_string(config, "mixture");
    int lv_setting = read_int(config, "lv_setting");

    YAML::Node data = config["data"];
    if (!data) {
        return entries;
    }

    auto build_entry = [&](const YAML::Node& node, const int layer) {
        MeasurementEntry entry;
        entry.name = read_string(node, "name");
        entry.measurement_type = measurement_type;
        entry.root_file = read_string(node, "root file");
        entry.mixture = mixture;
        entry.source = source;
        entry.filter = filter;
        entry.lv_setting = lv_setting;
        entry.hv = read_double(node, "hv");
        entry.other_hv = read_double(node, "other_hv");
        entry.scanned_hv = read_double(node, "scanned_hv");
        entry.layer = layer;
        if (!entry.name.empty()) {
            entries.push_back(std::move(entry));
        }
    };

    if (data.IsSequence()) {
        for (const auto& node : data) {
            build_entry(node, 0);
        }
    } else if (data.IsMap()) {
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

    return entries;
}
} // namespace ConfigUtils