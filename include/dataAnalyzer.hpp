#pragma once

#include <filesystem>
#include <iostream>

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include "types.hpp"
#include "hit.hpp"
#include "event.hpp"
#include "cluster.hpp"
#include "track.hpp"
#include "constants.hpp"


// ============================================================
// DataAnalyzer Class: Main data analysis section
// ============================================================
/// @class DataAnalyzer
/// @brief Main class for analyzing DCT data and producing summary statistics
class DataAnalyzer {
private:
    std::string config_path;
    std::filesystem::path output_directory = std::filesystem::path("data/output");
    std::vector<MeasurementEntry> parsed_entries;
    std::vector<SummaryStats> summaries;

public:
    explicit DataAnalyzer(const std::string& config_file_path);
    ~DataAnalyzer();

    /// @brief Main function to produce summary statistics for each measurement entry in the config file
    void produceSummaryStats();

    /// @brief Function to produce per-file relevant statistics for each measurement entry
    void producePerFileStats();

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
};