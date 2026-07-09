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
#include <functional>

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
#include "TLegend.h"

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

// ------------------------------------------------------------------------------------------
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
/// @namespace PlotterHelpers
/// @brief Collection of helper functions for plotting and styling ROOT objects
namespace PlotterHelpers {

    // --------------------------------------------------------------------------------------
    /// @namespace ATLASStyler
    /// @brief Collection of helper functions for applying ATLAS styling to plots
    namespace ATLASStyler {

        /// @brief Helper function to draw the ATLAS label on a plot
        /// @param x_offset The horizontal offset for the label position
        /// @param y_offset The vertical offset for the label position
        /// @param status The status text to display (e.g., "Work in Progress")
        void drawATLASLabel(double x_offset, double y_offset, const std::string& status);

        /// @brief Helper function to draw a title on a plot
        /// @param obj The ROOT object (e.g., TGraph, TH1) to which the title will be applied
        /// @param x_offset The horizontal offset for the title position
        /// @param y_offset The vertical offset for the title position
        void drawPlotTitle(TObject* obj, double x_offset, double y_offset);

        /// @brief Helper function to draw a legend on a plot
        /// @param mg A pointer to the TMultiGraph object for which the legend will be drawn
        /// @param x_offset The horizontal offset for the legend position
        /// @param y_offset The vertical offset for the legend position
        void drawATLASLegend(TMultiGraph* mg, double x_offset, double y_offset);

        /// @brief Helper function to adjust the color bar of a 2D histogram dynamically
        /// based on the maximum value in the histogram
        /// @param h2 A pointer to the TH2 histogram object
        /// @param pad A pointer to the TPad object on which the histogram is drawn
        void adjustDynamicCB(TH2* h2, TPad* pad);

        /// @brief Helper function to apply ATLAS styling to a given plot object and canvas
        /// @param obj The ROOT object (e.g., TGraph, TH1) to which the style will be applied
        /// @param pad A pointer to the TPad object on which the object is drawn (optional)
        void applyATLASStyle(TObject* obj, TPad* pad = nullptr);
    }   // namespace ATLASStyler

    // --------------------------------------------------------------------------------------
    /// @namespace PlotStyler
    /// @brief Namespace for functions that apply specific styling to different ROOT object types
    namespace PlotStyler {

        using PlotStylerFunc = std::function<void(TObject*, TCanvas*, TClass*)>;

        /// @brief Helper function to style efficiency vs HV plots
        /// @param obj The ROOT object (e.g., TGraph, TH1) to which the style will be applied
        /// @param canvas The TCanvas on which the object is drawn
        void styleEfficiencyVsHV(TObject* obj, TCanvas* canvas);

        /// @brief Default styling function for generic plots
        /// @param obj The ROOT object (e.g., TGraph, TH1) to which the style will be applied
        /// @param canvas The TCanvas on which the object is drawn
        /// @param class_type The TClass pointer representing the type of the ROOT object
        void styleDefaultPlot(TObject* obj, TCanvas* canvas, TClass* class_type);

        /// @brief Helper function to get a registry of plot styling functions for different
        /// ROOT object types
        /// @return A map of TClass pointers to corresponding styling functions
        std::unordered_map<PlotCategory, PlotStylerFunc> getPlotStyleRegistry();
    }   // namespace PlotStyler

    /// @brief Helper function to automatically export all relevant plots from a ROOT file
    /// to a specified directory in PDF format, applying ATLAS styling
    /// @param root_file_path The filesystem path to the input ROOT file
    /// @param target_plots_dir The filesystem path to the directory where the PDF files
    /// should be saved
    void autoExportToATLASPDF(
        const std::string& root_file_path,
        const std::filesystem::path& target_plots_dir
    );

    /// @brief Helper function to build global TMultiGraph objects for all metrics across
    /// all layers in a given configuration directory, and save them to the specified output path
    /// @param config_dir A pointer to the TDirectory representing the configuration directory 
    /// in the ROOT file
    /// @param config_output_path The filesystem path to the directory where the output files
    /// should be saved
    void buildGlobalMultiGraphs(
        TDirectory* config_dir,
        const std::filesystem::path& config_output_path
    );
}

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
