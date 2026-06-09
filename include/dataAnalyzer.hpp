#pragma once

#include <algorithm>
#include <filesystem>
#include <iostream>

#include <string>
#include <array>
#include <vector>
#include <unordered_set>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "types.hpp"
#include "hit.hpp"
#include "event.hpp"
#include "cluster.hpp"
#include "track.hpp"
#include "constants.hpp"


// ==========================================================================================
// Analysis utility functions for plotting and calculating statistics from processed DCT data
// ==========================================================================================
/// @namespace Utilities
/// @brief Namespace for utility functions used in the DataAnalyzer class
namespace Utilities {

    /// @struct ColumnsShift
    /// @brief Struct to define the start and end of a column and the shift to apply
    /// for remapping raw strip numbers to a continuous range for plotting
    struct ColumnsShift {
        int start;
        int end;
        int shift;
    };

    /// @brief Constant array defining the column shifts for remapping raw strip numbers
    /// to a continuous range for plotting
    constexpr std::array<ColumnsShift, 1> columnShifts = {{
        {16, 31, -8}    // Shifting strips 16-31 down by 8 to create a continuous range of 0-23 for plotting
    }};

    /// @brief Utility function to setup the directory structure in the input ROOT file
    /// @param input_file Pointer to the input ROOT file where the directory structure will be created
    void setupDirectory(TFile* input_file);

    /// @brief Utility function to remap raw strip numbers to a continuous range for plotting
    /// @param rawStrip The raw strip number to be remapped
    /// @return The remapped strip number for plotting
    int remapStrip(int rawStrip);
}

/// @namespace perFileHelpers
/// @brief Namespace for helper plotting functions producing per-filerelevant statistics
/// for each measurement entry
namespace perFileHelpers {
    void plotDtVsStrip(TFile* input_file, std::string entry_name);
}


// ==========================================================================================
// DataAnalyzer Class: Main data analysis section
// ==========================================================================================
/// @class DataAnalyzer
/// @brief Main class for analyzing DCT data and producing summary statistics
class DataAnalyzer {
public:
    explicit DataAnalyzer(const std::string& config_file_path);
    ~DataAnalyzer();

    /// @brief Main function to produce summary statistics for each measurement entry in the
    /// config file
    void produceSummaryStats();

    /// @brief Function to produce per-file relevant statistics for each measurement entry
    /// @param input_file Pointer to the input ROOT file containing processed DCT data for a
    /// specific measurement entry
    void producePerFileStats(TFile* input_file);

    // Getters
    const std::vector<MeasurementEntry>& getEntries() const { return parsed_entries; }
    const std::vector<SummaryStats>& getSummaries() const { return summaries; }

    // Utility funtions
    void applyCut();

    /* WIP: Move all this into a plotter class
    void plotEfficiencies();
    void plotRawDataDistributions();
    void plotStripDistributions();
    void plotDtDistributions();
    void plotDtvsStripDistributions();
    void plotToTDistributions();
    void plotToTvsStripDistributions();

    void plotClusterSizeDistributions();

    void plotTrackSizeDistributions();
    void plotTrackLengthDistributions();
    void plotTrackWidthDistributions();

    // Utility plotting functions
    void addATLASLogo(TCanvas* canvas);

    // Canvas plotting functions
    void makeCanvas(TGraph* graph, const std::string& title, const std::string& x_label, const std::string& y_label, const std::string& output_filename);
    */

private:
    std::string _config_path;
    std::filesystem::path _output_directory = std::filesystem::path("data/output");
    std::vector<MeasurementEntry> parsed_entries;
    std::vector<SummaryStats> summaries;
};