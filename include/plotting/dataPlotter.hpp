#pragma once

#include <string>
#include <vector>
#include <array>
#include <map>
#include <filesystem>

#include "core/types.hpp"

class TFile;
class TDirectory;

/// @namespace Utilities
/// @brief Namespace for general utility functions used across the plotting code
namespace Utilities {

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

    TFile* initializeAnalysisFile();
    TDirectory* setupScanDirectories(TDirectory* config_dir, const std::string& group_name);

    std::map<std::string, MetricsData> extractScanData(const std::string& summary_file_path);
    std::map<std::string, std::vector<FitResult>> extractCrossGroupFits(TDirectory* base_dir);

    void plotGlobalMetrics(TDirectory* scan_dir, const MetricsData& scan_data);
    void plotLayerMetrics(TDirectory* scan_dir, const std::map<std::string, LayerSeries>& layer_metrics);
    void plotStripMetrics(TDirectory* scan_dir, const std::map<std::string, std::map<int, std::map<int, StripSeries>>>& strip_metrics);
    void plotCrossGroupFits(TDirectory* config_dir, const std::map<std::string, std::vector<FitResult>>& all_fits);

    void cumulativeAnalysisRootFile();
    void cumulativeAnalysisPlots();

private:
    std::filesystem::path _output_directory;
    std::filesystem::path _analysis_root_file;

    /// @brief Map to hold parsed ConfigData structs for each scan, indexed by the configuration
    /// file path
    std::map<std::string, ConfigData> _parsed_configs;
};
