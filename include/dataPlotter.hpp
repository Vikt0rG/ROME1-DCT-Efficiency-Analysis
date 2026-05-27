#pragma once

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string_view>

#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>

#include <TGraph.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "types.hpp"
#include "configParser.hpp"

/// TODO: 
/// - Need function to style graphs (markers, colors, ticks, etc.) / ATLAS Style
/// - ATLAS Logo

// ============================================================
// DataPlotter Class: Plotting summary statistics
// ============================================================
class DataPlotter {
private:
    std::string config_path;
    std::filesystem::path output_directory = std::filesystem::path("data/output");
    std::vector<MeasurementEntry> parsed_entries;
    std::string summary_root_file;

public:
    explicit DataPlotter(const std::string& config_file_path);

    void produceSummaryPlots();

    // Accessors
    const TGraph* getGraphForMetric(const std::string& metric_name) const;

    // Utility functions
    static void removeLayerFromGraph(TGraph* graph, int layer_to_remove);
};
