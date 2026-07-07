#pragma once

#include <iostream>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string_view>

#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include "TClass.h"
#include "TKey.h"
#include <TGraph.h>
#include <TMultiGraph.h>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include "TLatex.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TSystem.h"

#include "types.hpp"
#include "utils.hpp"
#include "configParser.hpp"

// ==========================================================================================
// Plotting utility namespaces for creating graphs from the DataAnalyzer summary ROOT files
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
        std::array<std::vector<double>, 3> y_errors;
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
}

/// @enum PlotCategory
/// @brief Enum to represent different categories of plots that can be generated
enum class PlotCategory {
    Efficiency,
    EfficiencyVsHV,
    StripDistribution,
    DtVsStrip,
    ToT,
    ToTVsStrip,
    MultiplicityVsStrip,
    DelayVsStrip,
    ClusterSize,
    MeanClusterSize,
    Default
};

// ------------------------------------------------------------------------------------------
/// @brief Namespace for plotting helper functions used in the DataPlotter class
namespace PlotterHelpers {

    /// @brief Helper function to automatically export all relevant plots from a ROOT file
    /// to a specified directory in PDF format, applying ATLAS styling
    /// @param root_file_path The filesystem path to the input ROOT file
    /// @param target_plots_dir The filesystem path to the directory where the PDF files
    /// should be saved
    void autoExportToATLASPDF(
        const std::string& root_file_path,
        const std::filesystem::path& target_plots_dir
    );

    void removeLayerFromGraph(TGraph* graph, int layer_to_remove);

}

// ==========================================================================================
// DataPlotter Class: Plotting summary statistics
// ==========================================================================================
class DataPlotter {
public:
    DataPlotter(const std::vector<std::string>& config_paths, const std::filesystem::path& output_directory);

    void produceSummaryPlots();
    void exportPlotsToATLASPDF();

private:
    std::filesystem::path _output_directory;

    /// @brief Map to hold parsed ConfigData structs for each scan, indexed by the configuration
    /// file path
    std::map<std::string, ConfigData> _parsed_configs;

    /// @brief String path to the summary ROOT of a scan
    std::string _summary_root_file;
};
