#pragma once

#include <string>
#include <filesystem>

class TDirectory;

/// @namespace PlotterHelpers
/// @brief Namespace containing helper functions for plot styling and batch exporting plots
namespace PlotterHelpers {

/// @namespace BatchExporter
/// @brief Namespace containing helper functions for batch exporting plots from ROOT files
/// to PDF format with ATLAS styling
namespace BatchExporter {

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

}   // namespace BatchExporter
}   // namespace PlotterHelpers