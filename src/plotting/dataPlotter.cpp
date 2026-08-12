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
    return std::string(timestamp);
}
}   // namespace Utilities

// Anonymous namespace for metric names and other constants used in DataPlotter implementation
namespace {

    const std::vector<std::string> global_metrics = []() {
        std::vector<std::string> metrics;
        for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
            metrics.push_back("avg_time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta1");
            metrics.push_back("avg_time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta2");
            metrics.push_back("time_resolution_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta1");
            metrics.push_back("time_resolution_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta2");
        }
        return metrics;
    }();

    const std::vector<std::string> layer_metrics = {
        "eff_eta1_external", "eff_eta2_external", "eff_or_external", "eff_and_external",
        "eff_eta1_rpc", "eff_eta2_rpc", "eff_or_rpc", "eff_and_rpc",
        "track_eff_eta1_external", "track_eff_eta2_external", "track_eff_or_external", "track_eff_and_external",
        "track_eff_eta1_rpc", "track_eff_eta2_rpc", "track_eff_or_rpc", "track_eff_and_rpc",
        "avg_cluster_size_eta1_layers", "avg_cluster_size_eta2_layers",
        "noise_rate_eta1", "noise_rate_eta2"
    };

    const std::vector<std::string> scalar_metrics = {
        "avg_cluster_size_eta1", "avg_cluster_size_eta2", "noise_rate"
    };

    const std::vector<std::string> strip_layer_metrics = {
        "noise_rate_strips_eta1", "noise_rate_strips_eta2",
        "avg_tot_eta1", "avg_tot_eta2",
        "track_avg_tot_eta1", "track_avg_tot_eta2",
        "avg_multiplicity_eta1", "avg_multiplicity_eta2",
        "track_avg_multiplicity_eta1", "track_avg_multiplicity_eta2"
    };

    const std::vector<std::string> raw_metrics = []() {
        std::vector<std::string> metrics;
        for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
            metrics.push_back("time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta1");
            metrics.push_back("time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta2");
        }
        return metrics;
    }();

}   // anonymous namespace

// ==========================================================================================
// DataPlotter class implementation for plotting summary statistics
// ==========================================================================================
DataPlotter::DataPlotter(
    const std::vector<std::string>& config_paths,
    const std::filesystem::path& output_directory)
{
    _output_directory = output_directory;
    _parsed_configs = Utilities::parseConfigs(config_paths);
}

TFile* DataPlotter::initializeAnalysisFile() {
    std::filesystem::create_directories(_output_directory / "analysis");

    _analysis_root_file = _output_directory / "analysis" / (Utilities::getTimestamp() + "_analysis.root");
    TFile* analysis_root = TFile::Open(_analysis_root_file.c_str(), "RECREATE");
    PathUtils::verifyROOTFile(analysis_root, _analysis_root_file.string());

    std::vector<std::string> dir_names;
    for (const auto& [config_path, _] : _parsed_configs) {
        dir_names.push_back(std::filesystem::path(config_path).stem().string());
    }
    PathUtils::setupDirectories(analysis_root, dir_names);
    return analysis_root;
}

TDirectory* DataPlotter::setupScanDirectories(TDirectory* config_dir, const std::string& grp_name) {
    if (!config_dir) return nullptr;

    std::string scan_group_str = "group_" + grp_name;
    TDirectory* scan_dir = PathUtils::ensureDirectory(config_dir, scan_group_str.c_str());

    if (!scan_dir) return nullptr;

    PathUtils::ensureDirectory(scan_dir, "efficiency_analysis");
    PathUtils::ensureDirectory(scan_dir, "cluster_analysis");
    PathUtils::ensureDirectory(scan_dir, "noise_rate_layers_analysis");
    PathUtils::ensureDirectory(scan_dir, "noise_rate_strips_analysis");
    PathUtils::ensureDirectory(scan_dir, "tot_analysis");
    PathUtils::ensureDirectory(scan_dir, "multiplicity_analysis");
    PathUtils::ensureDirectory(scan_dir, "tof_analysis");

    return scan_dir;
}

// Helper function to extract scan data from a summary ROOT file and return a map of MetricsData
std::map<std::string, DataPlotter::MetricsData> DataPlotter::extractScanData(
    const std::string& summary_file_path)
{
    std::map<std::string, MetricsData> result;

    TFile* summary_root_file = TFile::Open(summary_file_path.c_str(), "READ");
    if (!summary_root_file || summary_root_file->IsZombie()) {
        if (summary_root_file) { summary_root_file->Close(); delete summary_root_file; }
        throw std::runtime_error("Failed to open summary root file: " + summary_file_path);
    }

    TTree* summary_tree = summary_root_file->Get<TTree>("summary");
    if (!summary_tree) {
        summary_root_file->Close(); delete summary_root_file;
        throw std::runtime_error("Summary tree not found in summary root file.");
    }

    TTreeReader readerSummary(summary_tree);

    TTreeReaderValue<std::string> group_name(readerSummary, "group_name");
    TTreeReaderValue<int> scanned_layer(readerSummary, "scanned_layer");
    TTreeReaderValue<double> scanned_hv(readerSummary, "scanned_hv");
    TTreeReaderValue<double> other_hv(readerSummary, "other_hv");

    // Readers for Scalar Metrics
    std::vector<std::unique_ptr<TTreeReaderValue<double>>> scalar_values;
    for (const auto& name : scalar_metrics) {
        scalar_values.push_back(std::make_unique<TTreeReaderValue<double>>(readerSummary, name.c_str()));
    }

    // Readers for Global Metrics
    std::vector<std::unique_ptr<TTreeReaderValue<double>>> global_values;
    std::vector<std::unique_ptr<TTreeReaderArray<double>>> global_errors;
    for (const auto& name : global_metrics) {
        global_values.push_back(std::make_unique<TTreeReaderValue<double>>(readerSummary, name.c_str()));
        global_errors.push_back(std::make_unique<TTreeReaderArray<double>>(readerSummary, (name + "_error").c_str()));
    }

    // Readers for Layer Metrics
    std::vector<std::unique_ptr<TTreeReaderArray<double>>> layer_arrays;
    std::vector<std::unique_ptr<TTreeReaderArray<double>>> layer_error_arrays;
    for (const auto& name : layer_metrics) {
        layer_arrays.push_back(std::make_unique<TTreeReaderArray<double>>(readerSummary, name.c_str()));
        layer_error_arrays.push_back(std::make_unique<TTreeReaderArray<double>>(readerSummary, (name + "_error").c_str()));
    }

    // Readers for Strip Metrics
    std::vector<std::unique_ptr<TTreeReaderArray<double>>> strip_arrays;
    std::vector<std::unique_ptr<TTreeReaderArray<double>>> strip_error_arrays;
    for (const auto& name : strip_layer_metrics) {
        strip_arrays.push_back(std::make_unique<TTreeReaderArray<double>>(readerSummary, name.c_str()));
        strip_error_arrays.push_back(std::make_unique<TTreeReaderArray<double>>(readerSummary, (name + "_error").c_str()));
    }

    // Readers for Raw Metrics
    std::vector<std::unique_ptr<TTreeReaderValue<std::vector<int>>>> raw_arrays;
    for (const auto& name : raw_metrics) {
        raw_arrays.push_back(std::make_unique<TTreeReaderValue<std::vector<int>>>(readerSummary, name.c_str()));
    }

    while (readerSummary.Next()) {
        std::string grp_name = *group_name;
        int scan_lyr = *scanned_layer;
        double scan_hv = *scanned_hv;
        double oth_hv = *other_hv;

        auto& current_scan = result[grp_name];

        // Extract Scalar Metrics (0D)
        for (size_t i = 0; i < scalar_metrics.size(); ++i) {
            const auto& metric_name = scalar_metrics[i];
            current_scan.scalar_x[metric_name].push_back(scan_hv);
            current_scan.scalar_y[metric_name].push_back(**scalar_values[i]);
        }

        // Extract Global Metrics (1D Graphs)
        for (size_t i = 0; i < global_metrics.size(); ++i) {
            const auto& metric_name = global_metrics[i];

            // Check if the branch actually exists in the file
            if (global_values[i]->GetSetupStatus() == 0) {
                auto& series = current_scan.global_metrics[metric_name];

                series.x.push_back(scan_hv);
                series.y.push_back(**global_values[i]);

                // Read the low and high bounds from the [2] array
                series.y_errors_low.push_back((*global_errors[i])[0]);
                series.y_errors_high.push_back((*global_errors[i])[1]);
            }
        }

        // Extract Raw Metrics (2D Heatmaps)
        for (size_t i = 0; i < raw_metrics.size(); ++i) {
            const auto& metric_name = raw_metrics[i];

            if (raw_arrays[i]->GetSetupStatus() == 0) {
                // Store the entire vector mapped to the current High Voltage
                current_scan.raw_tof_data[metric_name][scan_hv] = **raw_arrays[i];
            }
        }
 
        // Extract Layer Metrics (1D)
        for (size_t i = 0; i < layer_metrics.size(); ++i) {
            const auto& metric_name = layer_metrics[i];
            const auto& vals = *layer_arrays[i];
            const auto& errs = *layer_error_arrays[i];

            for (int layer = 0; layer < LAYER_COUNT; ++layer) {
                if (std::isnan(vals[layer])) continue;

                double x_value = (layer == scan_lyr) ? scan_hv : oth_hv;

                auto& layer_series = current_scan.layer_metrics[metric_name];
                layer_series.x[layer].push_back(x_value);
                layer_series.y[layer].push_back(vals[layer]);

                layer_series.y_errors_low[layer].push_back(errs[2 * layer]);
                layer_series.y_errors_high[layer].push_back(errs[2 * layer + 1]);
            }
        }

        // Extract Strip Metrics (2D)
        for (size_t i = 0; i < strip_layer_metrics.size(); ++i) {
            const auto& metric_name = strip_layer_metrics[i];
            const auto& vals = *strip_arrays[i];
            const auto& errs = *strip_error_arrays[i];

            for (int layer = 0; layer < LAYER_COUNT; ++layer) {
                double x_value = (layer == scan_lyr) ? scan_hv : oth_hv;
                for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {

                    int flat_idx = layer * STRIPS_PER_LAYER + strip;

                    if (std::isnan(vals[flat_idx])) continue;

                    auto& strip_series = current_scan.strip_metrics[metric_name][layer][strip];
                    strip_series.x.push_back(x_value);
                    strip_series.y.push_back(vals[flat_idx]);
                    strip_series.y_error_low.push_back(errs[flat_idx]);
                    strip_series.y_error_high.push_back(errs[flat_idx]);
                }
            }
        }
    }

    summary_root_file->Close();
    delete summary_root_file;
    return result;
}

void DataPlotter::plotGlobalMetrics(TDirectory* scan_dir, const MetricsData& scan_data) {
    if (!scan_dir) return;

    TDirectory* tof_dir = scan_dir->GetDirectory("tof_analysis");
    if (!tof_dir) return;
    tof_dir->cd();

    // Build 1D Graphs (Avg ToF & Time Resolution)
    for (const auto& [metric_name, series] : scan_data.global_metrics) {
        if (series.x.empty()) continue;

        TGraphAsymmErrors* graph = new TGraphAsymmErrors(
            series.x.size(), series.x.data(), series.y.data(),
            nullptr, nullptr, series.y_errors_low.data(), series.y_errors_high.data()
        );

        graph->SetName(metric_name.c_str());
        graph->SetTitle((metric_name + ";HV [V];Value").c_str());
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kAzure + 2);
        graph->SetLineColor(kAzure + 2);

        graph->Write("", TObject::kOverwrite);
        delete graph;
    }

    // Build 2D Heatmaps (Raw ToF vs HV)
    for (const auto& [eta_side, hv_data_map] : scan_data.raw_tof_data) {
        if (hv_data_map.empty()) continue;

        std::string hist_name = "h2d_" + eta_side;
        std::string hist_title = eta_side + ";High Voltage [V];Time of Flight [Ticks];Entries";

        TH2D* heatmap = new TH2D(hist_name.c_str(), hist_title.c_str(), 
                                 15, 4500, 6000,
                                 14, -7, 7);

        for (const auto& [hv, tof_vector] : hv_data_map) {
            for (int tof : tof_vector) {
                heatmap->Fill(hv, tof);
            }
        }

        heatmap->Write("", TObject::kOverwrite);
        delete heatmap;
    }
}

void DataPlotter::plotLayerMetrics(
    TDirectory* scan_dir, const std::map<std::string, Utilities::LayerSeries>& layer_metrics) {
    if (!scan_dir) return;

    // Grab the pre-created subdirectories
    TDirectory* eff_dir  = scan_dir->GetDirectory("efficiency_analysis");
    TDirectory* clus_dir = scan_dir->GetDirectory("cluster_analysis");
    TDirectory* nois_dir = scan_dir->GetDirectory("noise_rate_layers_analysis");

    for (const auto& [metric_name, series] : layer_metrics) {

        // Route to the correct directory based on the metric name
        TDirectory* metric_dir = scan_dir; // default fallback
        if (metric_name.rfind("eff_", 0) == 0 || metric_name.rfind("track_eff_", 0) == 0) metric_dir = eff_dir;
        else if (metric_name.rfind("avg_cluster", 0) == 0) metric_dir = clus_dir;
        else if (metric_name.rfind("noise_rate_eta", 0) == 0) metric_dir = nois_dir;

        if (!metric_dir) continue;
        metric_dir->cd();

        TMultiGraph* multi_graph = new TMultiGraph();
        multi_graph->SetName(metric_name.c_str());
        multi_graph->SetTitle((metric_name + ";HV;Value").c_str());

        for (int layer = 0; layer < LAYER_COUNT; ++layer) {
            if (series.x[layer].empty()) continue;

            TGraphAsymmErrors* layer_graph = new TGraphAsymmErrors(
                series.x[layer].size(), series.x[layer].data(), series.y[layer].data(),
                nullptr, nullptr, series.y_errors_low[layer].data(), series.y_errors_high[layer].data()
            );

            // Here set a beautiful identifier for the graph to be used as a label in the legend
            layer_graph->SetName(Form("%s_layer%d", metric_name.c_str(), layer));
            layer_graph->SetTitle(Form("Layer %d", layer));
            layer_graph->SetMarkerStyle(20 + layer);
            layer_graph->SetMarkerColor(1 + layer);
            layer_graph->SetLineColor(1 + layer);
            multi_graph->Add(layer_graph, "P");

            // Save individual layer graphs in subfolders
            std::string layer_folder = "layer" + std::to_string(layer);
            if (TDirectory* l_dir = PathUtils::ensureDirectory(metric_dir, layer_folder.c_str())) {
                l_dir->cd();
                layer_graph->Write("", TObject::kOverwrite);
            }
        }

        // Save the combined multigraph
        metric_dir->cd();
        multi_graph->Write("", TObject::kOverwrite);
        delete multi_graph;
    }
}

void DataPlotter::plotStripMetrics(
    TDirectory* scan_dir, const std::map<std::string, std::map<int,
    std::map<int, Utilities::StripSeries>>>& strip_metrics) {

    if (!scan_dir) return;

    TDirectory* nois_dir = scan_dir->GetDirectory("noise_rate_strips_analysis");
    TDirectory* tot_dir = scan_dir->GetDirectory("tot_analysis");
    TDirectory* mult_dir = scan_dir->GetDirectory("multiplicity_analysis");

    for (const auto& [metric_name, layer_map] : strip_metrics) {

        // Route to the correct parent analysis directory (e.g., tot_analysis)
        TDirectory* metric_dir = scan_dir;
        if (metric_name.find("noise_rate_strips_eta") != std::string::npos) metric_dir = nois_dir;
        else if (metric_name.find("tot") != std::string::npos) metric_dir = tot_dir;
        else if (metric_name.find("multiplicity") != std::string::npos) metric_dir = mult_dir;

        if (!metric_dir) continue;

        // Loop through each layer to create a dedicated TMultiGraph
        for (const auto& [layer, strip_map] : layer_map) {

            std::string layer_folder = "layer" + std::to_string(layer);
            TDirectory* l_dir = PathUtils::ensureDirectory(metric_dir, layer_folder.c_str());
            if (!l_dir) continue;

            // Create the MultiGraph strictly for this layer (holding 24 strips)
            TMultiGraph* layer_multi_graph = new TMultiGraph();
            std::string mg_name = metric_name + "_layer" + std::to_string(layer);
            layer_multi_graph->SetName(mg_name.c_str());
            layer_multi_graph->SetTitle((mg_name + ";HV;Value").c_str());

            // Populate it with all the strips for this layer
            for (const auto& [strip, data] : strip_map) {
                if (data.x.empty()) continue;

                TGraphAsymmErrors* g = new TGraphAsymmErrors(
                    data.x.size(), data.x.data(), data.y.data(),
                    nullptr, nullptr, data.y_error_low.data(), data.y_error_high.data()
                );

                g->SetName(Form("%s_strip%d", metric_name.c_str(), strip));
                g->SetTitle(Form("Strip %d", strip));
                g->SetMarkerStyle(20 + (strip % 4));
                g->SetMarkerColor(1 + (strip % 9));
                g->SetLineColor(1 + (strip % 9));

                layer_multi_graph->Add(g, "P");

                // Switch into the layer folder to drop off the single-strip graph
                l_dir->cd();
                g->Write("", TObject::kOverwrite);
            }

            // Switch back up to the analysis folder (e.g., tot_analysis) 
            // to save the assembled MultiGraph alongside the layer folders
            metric_dir->cd();
            layer_multi_graph->Write("", TObject::kOverwrite);
            delete layer_multi_graph;
        }
    }
}

void DataPlotter::cumulativeAnalysisRootFile() {
    // Handle file and basic directories
    TFile* analysis_root = initializeAnalysisFile();

    for (const auto& [config_path, config_data] : _parsed_configs) {
        std::string config_name = std::filesystem::path(config_path).stem().string();
        TDirectory* config_dir = analysis_root->GetDirectory(config_name.c_str());

        // Extract scan data from the summary ROOT file
        auto scan_data_map = extractScanData(config_data.summary_root_file);

        // Plot metrics for each scan layer
        for (const auto& [group_name, scan_data] : scan_data_map) {
            TDirectory* scan_dir = setupScanDirectories(config_dir, group_name);

            plotGlobalMetrics(scan_dir, scan_data);
            plotLayerMetrics(scan_dir, scan_data.layer_metrics);
            plotStripMetrics(scan_dir, scan_data.strip_metrics);
        }
    }

    analysis_root->Close();
    delete analysis_root;
}

void DataPlotter::cumulativeAnalysisPlots() {
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
    std::cout << "[ATLAS Export] Completed rendering all global metrics successfully." << std::endl;
}