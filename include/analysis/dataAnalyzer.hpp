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
// Analysis utility/helper namespace for plotting and calculating statistics
// ==========================================================================================
/// @namespace summaryHelpers
/// @brief Namespace for helper functions to calculate summary statistics for each
/// measurement entry
namespace summaryHelpers {

    /// @struct Accumulator
    /// @brief Struct to accumulate sums and counts for calculating averages and errors
    struct Accumulator {
        double sum{0.0};  ///< Sum of values accumulated
        size_t hits{0};   ///< Count of hits accumulated
    };

    void getEfficiency(TFile* input_file,
        EfficiencyResults& efficiency_results, EfficiencyResults& efficiency_results_tracks);
    void getClusterSize(TFile* input_file, ClusterSizeResults& cluster_size_results);
    void getRate(TFile* input_file, RateResults& rate_results);
    void getDeadStrips(DeadStrips& dead_strips);
    void getAverageToT(TFile* input_file, ToTResults& tot_results, bool in_valid_track_only,
        bool in_beam_only, const DeadStrips& dead_strips, const RoI& region_of_interest);
    void getAverageMultiplicity(TFile* input_file, MultiplicityResults& mult_results,
        bool in_valid_track_only, DeadStrips& dead_strips);
    void processToF(TFile* input_file, ToFResults& tof_results,
        TimeResolutionResults& time_resolution_results, DeadStrips& dead_strips);
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
    /// @param data Reference to the MeasurementData structure for the current measurement entry
    void producePerFileStats(TFile* input_file, MeasurementData& data);

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