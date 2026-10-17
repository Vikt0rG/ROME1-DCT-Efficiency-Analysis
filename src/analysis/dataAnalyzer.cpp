#include <iostream>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <array>

#include <TFile.h>
#include <TTree.h>

#include "configParser.hpp"
#include "plotBatchExporter.hpp"

#include "utils.hpp"
#include "constants.hpp"

#include "analysis/dataAnalyzer.hpp"
#include "analysis/perFileHelpers.hpp"
#include "analysis/summaryHelpers.hpp"


// ==========================================================================================
// DataAnalyzer class implementation for analyzing processed DCT data and plotting results
// ==========================================================================================
DataAnalyzer::DataAnalyzer(const std::string& config_file_path, const std::string& output_directory_path)
    : _config_path({config_file_path})
    , _output_directory(output_directory_path) {
}

DataAnalyzer::~DataAnalyzer() {
}

void DataAnalyzer::producePerFileStats(TFile* input_file, MeasurementData& data, MeasurementMetadata& metadata) {

    TDirectory* old_dir = input_file->GetDirectory("analysis");
    if (old_dir) {
        input_file->Delete("analysis;*");
        input_file->rmdir("analysis");
    }

    perFileHelpers::_metadata = metadata;

    // Set up output directories for per-file statistics and plots
    std::vector<std::string> dir_names = {
        "analysis/efficiency_map",
        "analysis/strip",
        "analysis/cluster_size",
        "analysis/tot",
        "analysis/dt_strip",
        "analysis/tot_strip",
        "analysis/multiplicity_strip",
        "analysis/delay_strip",
        "analysis/tof"
    };
    PathUtils::setupDirectories(input_file, dir_names);

    // Produce relevant plots for this file using helper functions
    perFileHelpers::plotEfficiencyMap(input_file, USE_EXTERNAL_TRIGGER);
    perFileHelpers::plotStrip(input_file);
    perFileHelpers::plotCS(input_file);
    perFileHelpers::plotDtVsStrip(input_file, data.region_of_interest);
    perFileHelpers::plotToT(input_file, data.region_of_interest);
    perFileHelpers::plotToTVsStrip(input_file);
    perFileHelpers::plotToTVsDtVsStrip(input_file, data.region_of_interest);
    perFileHelpers::plotMultiplicity(input_file);
    perFileHelpers::plotDelay(input_file);
    perFileHelpers::plotToFs(input_file);

    // Write analysis directory into the file
    input_file->Write("", TObject::kOverwrite);

    // Export per-file plots to PDF
    std::string root_file_path = input_file->GetName();
    std::string target_plots_dir = (_output_directory / "plots" / std::filesystem::path(input_file->GetName()).stem()).string();
    PlotterHelpers::BatchExporter::autoExportToATLASPDF(root_file_path, target_plots_dir);
}

void DataAnalyzer::produceGlobalStats(TFile* input_file, MeasurementData& data) {

    using namespace summaryHelpers;

    getEfficiency(input_file, data.efficiency_results, data.efficiency_results_tracks);
    getClusterSize(input_file, data.cluster_size_results);
    getRate(input_file, data.rate_results, false);
    getRate(input_file, data.rate_results_tracks, true);
    getDeadStrips(data.dead_strips);
    getAverageToT(input_file, data.tot_results, false, false, data.dead_strips, data.region_of_interest);
    getAverageToT(input_file, data.tot_results_tracks, false, true, data.dead_strips, data.region_of_interest);
    getAverageToT(input_file, data.tot_results_beam, true, false, data.dead_strips, data.region_of_interest);
    getAverageToT(input_file, data.tot_results_beam_track, true, true, data.dead_strips, data.region_of_interest);
    getAverageMultiplicity(input_file, data.multiplicity_results, false, data.dead_strips);
    getAverageMultiplicity(input_file, data.multiplicity_results_tracks, true, data.dead_strips);
    processToF(input_file, data.tof_results, data.time_resolution_results, data.dead_strips);
}

void DataAnalyzer::produceSummaryStats() {

    // Process each config file and build the list of measurement entries and summaries
    std::cout << _output_directory << std::endl;
    std::filesystem::create_directories(_output_directory / "root_summaries");

    // Prepare a measurement structure (measurement's metadata + data/statistics) and a root
    // file for this config
    ScanData scan;
    scan.config_path = _config_path;
    scan.metadata = ConfigUtils::parseMeasurementMetadata(_config_path);

    std::filesystem::path config_stem = std::filesystem::path(_config_path).stem();
    std::filesystem::path summary_root_path = _output_directory / "root_summaries" / (config_stem.string() + "_summary.root");

    // Create a ROOT file and tree to store summary statistics for this config
    TFile summary_root_file(summary_root_path.string().c_str(), "RECREATE");
    TTree* summary_tree = new TTree("summary", "summary");

    MeasurementMetadata metadata;
    MeasurementData data;

    // Set up branches for the summary tree
    summaryHelpers::setupBranches(summary_tree, metadata, data);

    // Process each measurement entry in this config file
    for (const auto& metadata_entry : scan.metadata) {
        if (metadata_entry.root_file.empty()) {
            std::cout << "Warning: Skipping entry '" << metadata_entry.name << "' due to missing ROOT file path." << std::endl;
            continue;
        }

        TFile* input_file = TFile::Open(metadata_entry.root_file.c_str(), "UPDATE");
        if (!input_file || input_file->IsZombie()) {
            std::cout << "Error: Failed to open ROOT file '" << metadata_entry.root_file << "' for entry '" << metadata_entry.name << "'. Skipping this entry." << std::endl;
            if (input_file) {
                std::cout << "Closing and deleting invalid ROOT file object for entry '" << metadata_entry.name << "'." << std::endl;
                input_file->Close();
                delete input_file;
            }
            continue;
        }

        data.clear();
        metadata = metadata_entry;

        // Calculate and fill per-file relevant statistics for this measurement entry and save into the input ROOT file
        producePerFileStats(input_file, data, metadata);

        // Calculate and store file statistics into the summary tree for this measurement entry
        produceGlobalStats(input_file, data);

        // Fill the summary tree
        summary_root_file.cd();
        summary_tree->Fill();
        scan.data.push_back(std::move(data));

        input_file->Close();
        delete input_file;
    }

    summary_root_file.cd();
    summary_tree->Write("", TObject::kOverwrite);
    summary_root_file.Close();
}