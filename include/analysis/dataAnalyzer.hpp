#pragma once

#include <regex>
#include <string>
#include <array>
#include <vector>
#include <filesystem>
#include "core/types.hpp"

class TFile;
class TH2;

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
    /// @param data Reference to the MeasurementData structure for the current measurement entry
    /// @param metadata Reference to the MeasurementMetadata structure for the current measurement entry
    void producePerFileStats(TFile* input_file, MeasurementData& data, MeasurementMetadata& metadata);

    /// @brief Function to produce global statistics for each measurement entry in the config file
    /// @param input_file Pointer to the input ROOT file containing processed DCT data for a
    /// specific measurement entry
    /// @param data Reference to the MeasurementData structure for the current measurement entry
    void produceGlobalStats(TFile* input_file, MeasurementData& data);

private:
    std::string _config_path;
    std::filesystem::path _output_directory;
    std::vector<MeasurementData> _measurement_data;
};