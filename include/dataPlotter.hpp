#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "types.hpp"

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
