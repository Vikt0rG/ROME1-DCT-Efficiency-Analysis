#pragma once

#include <string>
#include <vector>
#include <array>
#include <map>
#include <filesystem>

#include "core/types.hpp"

class TFile;
class TDirectory;

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

    /// @struct StripSeries
    /// @brief Struct to hold x and y, as well as y_errors data for a specific metric across
    /// different strips in a layer
    struct StripSeries {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> y_error_low;
        std::vector<double> y_error_high;
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

    struct MetricsData {
        std::map<std::string, Utilities::LayerSeries> layer_metrics;
        std::map<std::string, std::map<int, std::map<int, Utilities::StripSeries>>> strip_metrics;

        // For the scalar metrics
        std::map<std::string, std::vector<double>> scalar_x;
        std::map<std::string, std::vector<double>> scalar_y;
    };

    TFile* initializeAnalysisFile();
    TDirectory* setupScanDirectories(TDirectory* config_dir, const std::string& group_name);
    void plotLayerMetrics(TDirectory* scan_dir, const std::map<std::string, Utilities::LayerSeries>& layer_metrics);
    void plotStripMetrics(TDirectory* scan_dir, const std::map<std::string, std::map<int, std::map<int, Utilities::StripSeries>>>& strip_metrics);
    std::map<std::string, MetricsData> extractScanData(const std::string& summary_file_path);

    void cumulativeAnalysisRootFile();
    void cumulativeAnalysisPlots();

private:
    std::filesystem::path _output_directory;
    std::filesystem::path _analysis_root_file;

    /// @brief Map to hold parsed ConfigData structs for each scan, indexed by the configuration
    /// file path
    std::map<std::string, ConfigData> _parsed_configs;
};
