#pragma once

#include <string>
#include <array>
#include <vector>
#include <filesystem>

#include "types.hpp"

class TFile;

// ==========================================================================================
// Analysis utility/helper namespace for plotting and calculating statistics
// ==========================================================================================
/// @namespace perFileHelpers
/// @brief Namespace for helper plotting functions producing per-filerelevant statistics
/// for each measurement entry
namespace perFileHelpers {

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
        {16, 31, -8} // Shifting strips 16-31 down by 8 to create a continuous range of 0-23
    }};

    /// @brief Helper function to plot strip 1D distributions for all hits and valid tracks
    /// @param input_file Pointer to the input ROOT file containing processed DCT data for
    /// a specific measurement entry
    void plotStrip(TFile* input_file);

    /// @brief Helper function to plot dt vs strip 2D distributions for all hits and valid
    /// tracks
    /// @param input_file Pointer to the input ROOT file containing processed DCT data for
    /// a specific measurement entry
    void plotDtVsStrip(TFile* input_file);

    /// @brief Helper function to plot ToT vs strip 2D distributions for all hits and valid
    /// tracks
    /// @param input_file Pointer to the input ROOT file containing processed DCT data for
    /// a specific measurement entry
    void plotToTVsStrip(TFile* input_file);

    /// @brief Helper function to plot multiplicity/delay vs strip 2D distributions for all
    /// hits and valid tracks
    ///
    /// Multiplicity is defined as the number of hits in the same event and in the same strip
    ///
    /// Delay is defined as the time difference w.r.t. the first hit of the event in a strip
    /// @param input_file Pointer to the input ROOT file containing processed DCT data for
    /// a specific measurement entry
    void plotMultiplicityAndDelayVsStrip(TFile* input_file);
}


// ==========================================================================================
// DataAnalyzer Class: Main data analysis section
// ==========================================================================================
/// @class DataAnalyzer
/// @brief Main class for analyzing DCT data and producing summary statistics
class DataAnalyzer {
public:
    explicit DataAnalyzer(const std::string& config_file_path, const std::string& output_directory_path);
    ~DataAnalyzer();

    /// @brief Main function to produce summary statistics for each measurement entry in the
    /// config file
    void produceSummaryStats();

    /// @brief Function to produce per-file relevant statistics for each measurement entry
    /// @param input_file Pointer to the input ROOT file containing processed DCT data for a
    /// specific measurement entry
    void producePerFileStats(TFile* input_file);

private:
    std::string _config_path;
    std::filesystem::path _output_directory;
    std::vector<MeasurementData> _measurement_data;
};