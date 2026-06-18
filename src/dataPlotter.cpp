#include "dataPlotter.hpp"

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
}

// ------------------------------------------------------------------------------------------
namespace PlotterHelpers {

// Anonymous namespace for internal helpers & ATLAS styling
namespace {
    void drawATLASLabel(double x, double y, const std::string& status) {
        TLatex l;
        l.SetNDC();
        l.SetTextFont(72); // Bold-Italic Helvetica for "ATLAS"
        l.SetTextColor(kBlack);
        l.SetTextSize(0.04);
        l.DrawLatex(x, y, "ATLAS");

        if (!status.empty()) {
            TLatex s;
            s.SetNDC();
            s.SetTextFont(42); // Standard Helvetica for the status text
            s.SetTextColor(kBlack);
            s.SetTextSize(0.04);
            // Offset slightly to the right of "ATLAS"
            s.DrawLatex(x + 0.13, y, status.c_str()); 
        }
    }

    void applyATLASStyle(TObject* obj) {
        TAxis* xAxis = nullptr;
        TAxis* yAxis = nullptr;
        TAxis* zAxis = nullptr;

        // Safely extract axes depending on the object type
        if (auto h = dynamic_cast<TH1*>(obj)) {
            xAxis = h->GetXaxis();
            yAxis = h->GetYaxis();
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.2);
            h->SetLineColor(kBlack);
            h->SetLineWidth(2);
            h->SetStats(0);
        } else if (auto g = dynamic_cast<TGraph*>(obj)) {
            xAxis = g->GetXaxis();
            yAxis = g->GetYaxis();
            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.2);
            g->SetLineColor(kBlack);
            g->SetLineWidth(2);
            g->SetStats(0);
        } else if (auto h2 = dynamic_cast<TH2*>(obj)) {
            xAxis = h2->GetXaxis();
            yAxis = h2->GetYaxis();
            zAxis = h2->GetZaxis();
            h2->SetStats(0);
        }

        // Apply ATLAS standard font rules (Font 42 = Helvetica, Font 43 = Precise Pixel size)
        if (xAxis && yAxis) {
            xAxis->SetLabelFont(42);
            xAxis->SetTitleFont(42);
            xAxis->SetLabelSize(0.04);
            xAxis->SetTitleSize(0.05);
            xAxis->SetTitleOffset(1.1);

            yAxis->SetLabelFont(42);
            yAxis->SetTitleFont(42);
            yAxis->SetLabelSize(0.04);
            yAxis->SetTitleSize(0.05);
            yAxis->SetTitleOffset(1.3);
        }

        if (zAxis) {
            zAxis->SetLabelFont(42);
            zAxis->SetTitleFont(42);
            zAxis->SetLabelSize(0.04);
            zAxis->SetTitleSize(0.05);
            zAxis->SetTitleOffset(1.2);
        }
    }

    // Process a single directory recursively
    void scanDirectory(TDirectory* dir, const std::string& output_dir_path) {
        // Loop over all keys (objects saved) in the current directory
        TIter next(dir->GetListOfKeys());
        TKey* key = nullptr;

        while ((key = static_cast<TKey*>(next()))) {
            TClass* cl = TClass::GetClass(key->GetClassName());
            if (!cl) continue;

            // CASE A: Subdirectory found -> Dive in recursively
            if (cl->InheritsFrom(TDirectory::Class())) {
                TDirectory* subDir = dynamic_cast<TDirectory*>(key->ReadObj());
                if (subDir) {
                    scanDirectory(subDir, output_dir_path);
                }
                continue;
            }

            // CASE B: Target visual object found
            if (cl->InheritsFrom(TH1::Class()) || cl->InheritsFrom(TGraph::Class())) {
                TObject* obj = key->ReadObj();
                if (!obj) continue;

                // Create a clean, standard ATLAS aspect-ratio canvas (600x600 or 800x600)
                TCanvas* canvas = new TCanvas("c_temp", "", 800, 600);
                canvas->SetLeftMargin(0.16);
                canvas->SetRightMargin(0.05);
                canvas->SetTopMargin(0.06);
                canvas->SetBottomMargin(0.14);

                // Handle 2D histograms differently (usually require colz palette)
                if (cl->InheritsFrom(TH2::Class())) {
                    canvas->SetRightMargin(0.14); // Make room for color palette
                    obj->Draw("COLZ");
                } else if (cl->InheritsFrom(TGraph::Class())) {
                    // Force axis generation for TGraphs using "AP" (Axis + Points)
                    obj->Draw("AP"); 
                } else {
                    // 1D Histograms draw with errors by default
                    obj->Draw("E1"); 
                }

                // Apply cosmetics & branding
                applyATLASStyle(obj);
                drawATLASLabel(0.19, 0.88, "Internal");

                // Construct clean output file string
                std::string safe_name = obj->GetName();
                std::string export_path = output_dir_path + "/" + safe_name + ".pdf";

                // Save and drop memory overhead
                canvas->SaveAs(export_path.c_str());
                
                delete canvas;
                delete obj;
            }
        }
    }
}   // Anonymous namespace

void autoExportToATLASPDF(const std::string& root_file_path, const std::string& target_plots_dir) {
    // Ensure the output folder exists
    gSystem->mkdir(target_plots_dir.c_str(), kTRUE);

    // Open file in strict READ mode
    TFile* file = TFile::Open(root_file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        if (file) delete file;
        throw std::runtime_error("Export pipeline failed to open input file: " + root_file_path);
    }

    std::cout << "[ATLAS Export] Beginning recursive sweep of " << root_file_path << "..." << std::endl;
    
    // Launch recursive extraction sweep starting from file root
    scanDirectory(file, target_plots_dir);

    std::cout << "[ATLAS Export] Completed cleanly. Plots saved to: " << target_plots_dir << std::endl;

    file->Close();
    delete file;
}

void removeLayerFromGraph(TGraph* graph, int layer_to_remove) {
    if (!graph) { 
        std::cerr << "Graph pointer is null." << std::endl;
        return; 
    }
    TMultiGraph* multi_graph = dynamic_cast<TMultiGraph*>(graph);
    if (!multi_graph) { 
        std::cerr << "Failed to cast graph to TMultiGraph." << std::endl;
        return; 
    }
    TGraph* layer_graph = nullptr;
    auto* list = multi_graph->GetListOfGraphs();
    if (!list) { 
        std::cerr << "MultiGraph has no list of graphs." << std::endl;
        return; 
    }
    for (int i = 0; i < list->GetSize(); ++i) {
        TGraph* g = dynamic_cast<TGraph*>(list->At(i));
        if (g && g->GetName() && std::string(g->GetName()).find("layer" + std::to_string(layer_to_remove)) != std::string::npos) {
            layer_graph = g;
            break;
        }
    }
    if (layer_graph) {
        list->Remove(layer_graph);
        delete layer_graph;
    } else {
        std::cerr << "Layer graph not found for layer: " << layer_to_remove << std::endl;
    }
}
}

// ==========================================================================================
// DataPlotter class implementation for plotting summary statistics
// ==========================================================================================
DataPlotter::DataPlotter(const std::vector<std::string>& config_paths)
    : _config_paths(config_paths) {
    _parsed_configs = Utilities::parseConfigs(config_paths);
}

// Produce summary plots based on the parsed entries and summary root file
void DataPlotter::produceSummaryPlots() {
    std::filesystem::create_directories(_output_directory);

    // Create an output path following "DD-MM-YYYY_HH-MM-SS_analysis.root" format
    std::filesystem::path analysis_root_path = _output_directory / (Utilities::getTimestamp() + "_analysis.root");
    TFile* analysis_root = TFile::Open(analysis_root_path.c_str(), "RECREATE");
    PathUtils::verifyROOTFile(analysis_root, analysis_root_path.string());

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
        "eff_eta1_external_error",
        "eff_eta2_external_error",
        "eff_or_external_error",
        "eff_and_external_error",
        "eff_eta1_rpc",
        "eff_eta2_rpc",
        "eff_or_rpc",
        "eff_and_rpc",
        "eff_eta1_rpc_error",
        "eff_eta2_rpc_error",
        "eff_or_rpc_error",
        "eff_and_rpc_error",
        "track_eff_eta1_external",
        "track_eff_eta2_external",
        "track_eff_or_external",
        "track_eff_and_external",
        "track_eff_eta1_external_error",
        "track_eff_eta2_external_error",
        "track_eff_or_external_error",
        "track_eff_and_external_error",
        "track_eff_eta1_rpc",
        "track_eff_eta2_rpc",
        "track_eff_or_rpc",
        "track_eff_and_rpc",
        "track_eff_eta1_rpc_error",
        "track_eff_eta2_rpc_error",
        "track_eff_or_rpc_error",
        "track_eff_and_rpc_error",

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
        layer_arrays.reserve(layer_metrics.size());
        for (const auto& name : layer_metrics) {
            layer_arrays.push_back(std::make_unique<TTreeReaderArray<double>>(readerSummary, name.c_str()));
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
                if (layer_series.find(metric_name) == layer_series.end()) {
                    layer_series.emplace(metric_name, Utilities::LayerSeries{});
                }
                for (int layer = 0; layer < 3; ++layer) {
                    const double x_value = (layer == scanned_layer_value) ? scanned_hv_value : other_hv_value;
                    layer_series[metric_name].x[layer].push_back(x_value);
                    layer_series[metric_name].y[layer].push_back(values[layer]);
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
                // Give it a highly unique name combining scan environment and metric name
                std::string mg_unique_name = "scanned_layer" + std::to_string(scan_layer) + "_" + metric_name;
                multi_graph->SetName(mg_unique_name.c_str());
                std::string title_str = metric_name + "_scanned_layer" + std::to_string(scan_layer) + ";HV;Value";
                multi_graph->SetTitle(title_str.c_str());

                // Inner Loop: Iterate through layers (0, 1, 2)
                for (int layer = 0; layer < 3; ++layer) {
                    const auto& xs = series.x[layer];
                    const auto& ys = series.y[layer];
                    if (xs.empty()) { continue; }
                    
                    TGraph* layer_graph = new TGraph(static_cast<int>(xs.size()), xs.data(), ys.data());
                    
                    // Name reflects the physical layer being recorded
                    std::string layer_graph_name = mg_unique_name + "_layer" + std::to_string(layer);
                    layer_graph->SetName(layer_graph_name.c_str());
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
}

/*
// Accessors
const TGraph* DataPlotter::getGraphForMetric(const std::string& metric_name) const {
    std::filesystem::path config_stem = std::filesystem::path(config_path).stem();
    std::filesystem::path analysis_root_path = _output_directory / (config_stem.string() + "_analysis.root");
    TFile* file = TFile::Open(analysis_root_path.string().c_str(), "READ");
    if (!file || file->IsZombie()) {
        if (file) {
            file->Close();
            delete file;
            file = nullptr;
        }
        std::cerr << "Failed to open analysis root file: " << analysis_root_path << std::endl;
        return nullptr;
    }
    TGraph* graph = dynamic_cast<TGraph*>(file->Get(metric_name.c_str()));
    if (!graph) {
        file->Close();
        delete file;
        std::cerr << "Graph not found for metric: " << metric_name << std::endl;
        return nullptr;
    }
    return graph;
}
*/