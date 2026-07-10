#include <iostream>
#include <unordered_map>
#include <stdexcept>

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
#include "TAxis.h"

#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

#include "utils.hpp"
#include "configParser.hpp"
#include "plotting/plotStyler.hpp"
#include "plotting/plotBatchExporter.hpp"
#include "plotting/dataPlotter.hpp"

// ==========================================================================================
// Plotting utility functions for creating summaries from the DataAnalyzer summary ROOT files
// ==========================================================================================
namespace Utilities {
std::map<std::string, ConfigData> parseConfigs(const std::vector<std::string>& config_paths) {
    std::map<std::string, ConfigData> entries_by_scan;
    
    for (const std::string& config_path : config_paths) {
        std::string summary_path;
        auto metadata = ConfigUtils::parseMeasurementMetadata(config_path, &summary_path);
        
        if (!metadata.empty()) {
            entries_by_scan[config_path] = ConfigData{std::move(metadata), summary_path};
        }
    }
    return entries_by_scan;
}

std::string getTimestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm;
    localtime_r(&now_time_t, &now_tm);
    char timestamp[20];
    std::strftime(timestamp, sizeof(timestamp), "%d-%m-%Y_%H-%M-%S", &now_tm);
    return std::string(timestamp);}
}   // namespace Utilities

// ==========================================================================================
// DataPlotter class implementation for plotting summary statistics
// ==========================================================================================
DataPlotter::DataPlotter(const std::vector<std::string>& config_paths, const std::filesystem::path& output_directory) {
    _output_directory = output_directory;
    _parsed_configs = Utilities::parseConfigs(config_paths);
}

// Produce summary plots based on the parsed entries and summary root file
void DataPlotter::produceSummaryPlots() {
    std::filesystem::create_directories(_output_directory / "analysis");

    // Create an output path following "DD-MM-YYYY_HH-MM-SS_analysis.root" format
    _analysis_root_file = _output_directory / "analysis" / (Utilities::getTimestamp() + "_analysis.root");
    TFile* analysis_root = TFile::Open(_analysis_root_file.c_str(), "RECREATE");
    PathUtils::verifyROOTFile(analysis_root, _analysis_root_file.string());

    // Set up analysis ROOT file
    std::vector<std::string> dir_names;
    for (const auto& [config_path, _] : _parsed_configs) {
        std::string config_stem = std::filesystem::path(config_path).stem().string();
        dir_names.push_back(config_stem);
    }
    PathUtils::setupDirectories(analysis_root, dir_names);

    // Prepare pre-defined metric lists to read from the summary files
    // Layer-specific pre-defined arrays of length 3 (one per layer)
    const std::vector<std::string> layer_metrics = {
        "eff_eta1_external",
        "eff_eta2_external",
        "eff_or_external",
        "eff_and_external",
        "eff_eta1_rpc",
        "eff_eta2_rpc",
        "eff_or_rpc",
        "eff_and_rpc",
        "track_eff_eta1_external",
        "track_eff_eta2_external",
        "track_eff_or_external",
        "track_eff_and_external",
        "track_eff_eta1_rpc",
        "track_eff_eta2_rpc",
        "track_eff_or_rpc",
        "track_eff_and_rpc",

        "avg_cluster_size_eta1_layers",
        "avg_cluster_size_eta2_layers",

        "noise_rate_eta1",
        "noise_rate_eta2"
    };
    // Detector-wide scalar metrics
    const std::vector<std::string> scalar_metrics = {
        "avg_cluster_size_eta1",
        "avg_cluster_size_eta2",
        "noise_rate"
    };

    // Iterate over the configs and their data to produce summary plots
    for (const auto& [config_path, config_data] : _parsed_configs) {

        // Switch to the directory for this config in the analysis ROOT file
        TDirectory* config_dir = analysis_root->GetDirectory(std::filesystem::path(config_path).stem().string().c_str());
        config_dir->cd();

        // Extract information from the ConfigData struct for this config
        _summary_root_file = config_data.summary_root_file;

        // Read summary ROOT file from the map for each config
        TFile* summary_root_file = TFile::Open(_summary_root_file.c_str(), "READ");
        if (!summary_root_file || summary_root_file->IsZombie()) {
            if (summary_root_file) {
                summary_root_file->Close();
                delete summary_root_file;
            }
            throw std::runtime_error("Failed to open summary root file: " + _summary_root_file);
        }
        TTree* summary_tree = dynamic_cast<TTree*>(summary_root_file->Get("summary"));
        if (!summary_tree) {
            summary_root_file->Close();
            delete summary_root_file;
            throw std::runtime_error("Summary tree not found in summary root file.");
        }

        // Prepare data structures to hold the series for each metric and layer
        std::map<int, std::map<std::string, Utilities::LayerSeries>> layer_series_by_scan;
        std::map<int, std::map<std::string, std::vector<double>>> scalar_x_by_scan;
        std::map<int, std::map<std::string, std::vector<double>>> scalar_y_by_scan;

        // Prepare TTreeReader and readers for the relevant branches
        TTreeReader readerSummary(summary_tree);
        TTreeReaderValue<std::string> name(readerSummary, "name");
        TTreeReaderValue<int> scanned_layer(readerSummary, "scanned_layer");
        TTreeReaderValue<double> scanned_hv(readerSummary, "scanned_hv");
        TTreeReaderValue<double> other_hv(readerSummary, "other_hv");

        std::vector<std::unique_ptr<TTreeReaderArray<double>>> layer_arrays;
        std::vector<std::unique_ptr<TTreeReaderArray<double>>> layer_error_arrays;
        layer_arrays.reserve(layer_metrics.size());
        layer_error_arrays.reserve(layer_metrics.size());
        for (const auto& name : layer_metrics) {
            layer_arrays.push_back(std::make_unique<TTreeReaderArray<double>>(readerSummary, name.c_str()));
            layer_error_arrays.push_back(std::make_unique<TTreeReaderArray<double>>(readerSummary, (name + "_error").c_str()));
        }

        std::vector<std::unique_ptr<TTreeReaderValue<double>>> scalar_values;
        scalar_values.reserve(scalar_metrics.size());
        for (const auto& name : scalar_metrics) {
            scalar_values.push_back(std::make_unique<TTreeReaderValue<double>>(readerSummary, name.c_str()));
        }

        // Iterate over the summary tree entries and populate the series data structures
        // x := HV, y := metric value (e.g. efficiency), separate series for each layer (0, 1, 2)
        while (readerSummary.Next()) {
            const double scanned_hv_value = *scanned_hv;
            const int scanned_layer_value = *scanned_layer;
            double other_hv_value = *other_hv;

            // Populate layer-specific series
            auto& layer_series = layer_series_by_scan[scanned_layer_value];
            for (size_t i = 0; i < layer_metrics.size(); ++i) {
                const auto& metric_name = layer_metrics[i];
                const auto& values = *layer_arrays[i];
                const auto& errors = *layer_error_arrays[i];

                if (layer_series.find(metric_name) == layer_series.end()) {
                    layer_series.emplace(metric_name, Utilities::LayerSeries{});
                }
                for (int layer = 0; layer < 3; ++layer) {
                    const double x_value = (layer == scanned_layer_value) ? scanned_hv_value : other_hv_value;
                    layer_series[metric_name].x[layer].push_back(x_value);
                    layer_series[metric_name].y[layer].push_back(values[layer]);
                    layer_series[metric_name].y_errors[layer].push_back(errors[layer]);
                }
            }

            // Populate scalar metrics
            auto& scalar_x = scalar_x_by_scan[scanned_layer_value];
            auto& scalar_y = scalar_y_by_scan[scanned_layer_value];
            for (size_t i = 0; i < scalar_metrics.size(); ++i) {
                const auto& metric_name = scalar_metrics[i];
                if (scalar_x.find(metric_name) == scalar_x.end()) {
                    scalar_x.emplace(metric_name, std::vector<double>{});
                    scalar_y.emplace(metric_name, std::vector<double>{});
                }
                scalar_x[metric_name].push_back(scanned_hv_value);
                scalar_y[metric_name].push_back(**scalar_values[i]);
            }
        }

        // Once we have all the series data organized by scanned layer and metric, we can create
        // TMultiGraphs for each metric and populate them with TGraphs for each layer
        // Outer Loop: Iterate over the "Scanned Layers" (e.g., scan_layer = "scanned_layer_0")
        for (const auto& [scan_layer, layer_series] : layer_series_by_scan) {

            // Dynamically create/get the directory for this specific High Voltage scan environment
            std::string scan_layer_str = "scanned_layer_" + std::to_string(scan_layer);
            TDirectory* scan_dir = PathUtils::ensureDirectory(config_dir, scan_layer_str.c_str());
            if (!scan_dir) continue;

            // Create the specialized analysis metric folders inside this scanned layer's directory
            TDirectory* efficiencies_dir = PathUtils::ensureDirectory(scan_dir, "efficiency_analysis");
            TDirectory* clusters_dir     = PathUtils::ensureDirectory(scan_dir, "cluster_analysis");
            TDirectory* noise_dir        = PathUtils::ensureDirectory(scan_dir, "noise_rate_analysis");

            // Helper to return the correct metric parent directory
            auto pick_metric_directory = [&](std::string_view metric_name) -> TDirectory* {
                if (metric_name.rfind("eff_", 0) == 0 || metric_name.rfind("track_eff_", 0) == 0) {
                    return efficiencies_dir ? efficiencies_dir : scan_dir;
                }
                if (metric_name.rfind("avg_cluster_size", 0) == 0) {
                    return clusters_dir ? clusters_dir : scan_dir;
                }
                if (metric_name.rfind("noise_rate", 0) == 0) {
                    return noise_dir ? noise_dir : scan_dir;
                }
                return scan_dir;
            };

            // Middle Loop: Iterate over metrics
            for (const auto& [metric_name, series] : layer_series) {
                TDirectory* metric_dir = pick_metric_directory(metric_name);
                if (!metric_dir) continue;

                // Change ROOT context to the metric directory to host the TMultiGraph
                metric_dir->cd();
                
                TMultiGraph* multi_graph = new TMultiGraph();
                multi_graph->SetName(metric_name.c_str());
                std::string title_str = metric_name + ";HV;Value";
                multi_graph->SetTitle(title_str.c_str());

                // Inner Loop: Iterate through layers (0, 1, 2)
                for (int layer = 0; layer < 3; ++layer) {
                    const auto& xs = series.x[layer];
                    const auto& ys = series.y[layer];
                    const auto& y_errors = series.y_errors[layer];
                    if (xs.empty()) { continue; }
                    
                    TGraphErrors* layer_graph = new TGraphErrors(
                        static_cast<int>(xs.size()),
                        xs.data(), 
                        ys.data(),
                        nullptr,
                        y_errors.data()
                    );
                    
                    // Name reflects the physical layer being recorded
                    layer_graph->SetName(metric_name.c_str());
                    layer_graph->SetMarkerStyle(20 + layer);
                    layer_graph->SetMarkerColor(1 + layer);
                    layer_graph->SetLineColor(1 + layer);
                    
                    // Add this line graph into our Scanned Layer's MultiGraph combo pack
                    multi_graph->Add(layer_graph, "P");

                    // Save the raw individual component graph into its own measurement layer subfolder
                    std::string layer_folder_name = "layer" + std::to_string(layer);
                    TDirectory* layer_sub_dir = PathUtils::ensureDirectory(metric_dir, layer_folder_name.c_str());
                    if (layer_sub_dir) {
                        layer_sub_dir->cd();
                        layer_graph->Write("", TObject::kOverwrite); 
                    }
                }

                // Wrap up: Write the completed TMultiGraph to the metric directory level
                metric_dir->cd();
                multi_graph->Write("", TObject::kOverwrite);
                
                // Memory cleanup for this loop cycle
                delete multi_graph; 
            }
        }
    }   // End of loop over configs
    analysis_root->Close();
    delete analysis_root;
}

void DataPlotter::exportPlotsToATLASPDF() {
    std::string timestamp = Utilities::getTimestamp();
    std::filesystem::path target_plots_dir = _output_directory / "plots";
    std::filesystem::create_directories(target_plots_dir);

    TFile* file = TFile::Open(_analysis_root_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Could not open analysis file for PDF rendering: " << _analysis_root_file << std::endl;
        if (file) delete file;
        return;
    }

    std::cout << "[ATLAS Export] Automatically fetching and exporting all metrics..." << std::endl;

    // Loop through all top-level configuration directories
    TIter next_config(file->GetListOfKeys());
    TKey* config_key = nullptr;

    while ((config_key = static_cast<TKey*>(next_config()))) {
        TClass* cl = TClass::GetClass(config_key->GetClassName());
        if (!cl || !cl->InheritsFrom(TDirectory::Class())) continue;

        TDirectory* config_dir = dynamic_cast<TDirectory*>(config_key->ReadObj());
        if (!config_dir) continue;

        std::string config_name = config_key->GetName();
        std::filesystem::path config_output_path = target_plots_dir / config_name;

        std::cout << "  -> Processing configuration: " << config_name << std::endl;

        // Dynamically process every analysis directory and metric present in the file
        PlotterHelpers::BatchExporter::buildGlobalMultiGraphs(config_dir, config_output_path);

        delete config_dir;
    }

    file->Close();
    delete file;
    std::cout << "[ATLAS Export] Completed rendering all intrinsic metrics successfully." << std::endl;
}