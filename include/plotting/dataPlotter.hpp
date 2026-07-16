#pragma once

#include <string>
#include <vector>
#include <array>
#include <map>
#include <filesystem>

#include "core/types.hpp"

// ==========================================================================================
// Utilities Namespace: General utility functions for plotting
// ==========================================================================================
/// @namespace Utilities
/// @brief Namespace for general utility functions used across the plotting code
namespace Utilities {
    /// @struct LayerSeries
    /// @brief Struct to hold x and y, as well as y_errors data for a specific metric across
    /// different three layers
    struct LayerSeries {
        std::array<std::vector<double>, 3> x;
        std::array<std::vector<double>, 3> y;
        std::array<std::vector<double>, 3> y_errors_low;
        std::array<std::vector<double>, 3> y_errors_high;
    };
    
    /// @brief Utility function to parse measurement entries from provided configuration paths
    /// @param config_paths A vector of strings representing paths to YAML configuration files
    /// for different measurement entries
    /// @return A map of scan data indexed by the configuration file path
    std::map<std::string, ConfigData> parseConfigs(
        const std::vector<std::string>& config_paths
    );

    /// @brief Utility function to get a timestamp string for naming output files
    /// @return A string with the current timestamp in the format "DD-MM-YYYY_HH-MM-SS"
    std::string getTimestamp();
}   // namespace Utilities

// ==========================================================================================
// DataPlotter Class: Plotting summary statistics
// ==========================================================================================
class DataPlotter {
public:
    DataPlotter(
        const std::vector<std::string>& config_paths,
        const std::filesystem::path& output_directory
    );

    void produceSummaryPlots();
    void exportPlotsToATLASPDF();

private:
    std::filesystem::path _output_directory;
    std::filesystem::path _analysis_root_file;

    /// @brief Map to hold parsed ConfigData structs for each scan, indexed by the configuration
    /// file path
    std::map<std::string, ConfigData> _parsed_configs;

    /// @brief String path to the summary ROOT of a scan
    std::string _summary_root_file;
};
