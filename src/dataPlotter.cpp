#include "dataPlotter.hpp"

// Constructor implementation: Parse config file and populate parsed_entries vector
DataPlotter::DataPlotter(const std::string& config_file_path)
    : config_path(config_file_path) {
    parsed_entries = ConfigUtils::parseMeasurementEntries(config_file_path, &summary_root_file);
}

// Produce summary plots based on the parsed entries and summary root file
void DataPlotter::produceSummaryPlots() {
    std::filesystem::create_directories(output_directory);

    std::filesystem::path config_stem = std::filesystem::path(config_path).stem();
    std::filesystem::path analysis_root_path = output_directory / (config_stem.string() + "_analysis.root");

    // Set up analysis ROOT file
    TFile* summary_root = TFile::Open(summary_root_file.c_str(), "READ");
    if (!summary_root || summary_root->IsZombie()) {
        if (summary_root) {
            summary_root->Close();
            delete summary_root;
        }
        throw std::runtime_error("Failed to open summary root file: " + summary_root_file);
    }

    TTree* summary_tree = dynamic_cast<TTree*>(summary_root->Get("summary"));
    if (!summary_tree) {
        summary_root->Close();
        delete summary_root;
        throw std::runtime_error("Summary tree not found in summary root file.");
    }

    struct LayerSeries {
        std::array<std::vector<double>, 3> x;
        std::array<std::vector<double>, 3> y;
    };

    // Those are pre-defined arrays of length 3 (for 3 layers)
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

    // And those are scalar metrics that are detector-wide (not layer-specific)
    const std::vector<std::string> scalar_metrics = {
        "avg_cluster_size_eta1",
        "avg_cluster_size_eta2",
        "noise_rate"
    };

    std::map<int, std::map<std::string, LayerSeries>> layer_series_by_scan;
    std::map<int, std::map<std::string, std::vector<double>>> scalar_x_by_scan;
    std::map<int, std::map<std::string, std::vector<double>>> scalar_y_by_scan;

    std::unordered_map<std::string, double> other_hv_by_name;
    other_hv_by_name.reserve(parsed_entries.size());
    for (const auto& entry : parsed_entries) {
        other_hv_by_name.emplace(entry.name, entry.other_hv);
    }

    TTreeReader reader(summary_tree);
    TTreeReaderValue<std::string> entry_name(reader, "name");
    TTreeReaderValue<int> entry_layer(reader, "layer");
    TTreeReaderValue<double> hv(reader, "hv");

    std::vector<std::unique_ptr<TTreeReaderArray<double>>> layer_arrays;
    layer_arrays.reserve(layer_metrics.size());
    for (const auto& name : layer_metrics) {
        layer_arrays.push_back(std::make_unique<TTreeReaderArray<double>>(reader, name.c_str()));
    }

    std::vector<std::unique_ptr<TTreeReaderValue<double>>> scalar_values;
    scalar_values.reserve(scalar_metrics.size());
    for (const auto& name : scalar_metrics) {
        scalar_values.push_back(std::make_unique<TTreeReaderValue<double>>(reader, name.c_str()));
    }

    while (reader.Next()) {
        const double hv_value = *hv;
        double other_hv_value = hv_value;
        auto other_it = other_hv_by_name.find(*entry_name);
        if (other_it != other_hv_by_name.end() && other_it->second > 0.0) {
            other_hv_value = other_it->second;
        }
        const int scanned_layer = *entry_layer;
        auto& layer_series = layer_series_by_scan[scanned_layer];
        auto& scalar_x = scalar_x_by_scan[scanned_layer];
        auto& scalar_y = scalar_y_by_scan[scanned_layer];
        for (size_t i = 0; i < layer_metrics.size(); ++i) {
            const auto& metric_name = layer_metrics[i];
            const auto& values = *layer_arrays[i];
            if (layer_series.find(metric_name) == layer_series.end()) {
                layer_series.emplace(metric_name, LayerSeries{});
            }
            for (int layer = 0; layer < 3; ++layer) {
                const double x_value = (layer == scanned_layer) ? hv_value : other_hv_value;
                layer_series[metric_name].x[layer].push_back(x_value);
                layer_series[metric_name].y[layer].push_back(values[layer]);
            }
        }
        for (size_t i = 0; i < scalar_metrics.size(); ++i) {
            const auto& metric_name = scalar_metrics[i];
            if (scalar_x.find(metric_name) == scalar_x.end()) {
                scalar_x.emplace(metric_name, std::vector<double>{});
                scalar_y.emplace(metric_name, std::vector<double>{});
            }
            scalar_x[metric_name].push_back(hv_value);
            scalar_y[metric_name].push_back(**scalar_values[i]);
        }
    }

    TFile output_root(analysis_root_path.string().c_str(), "RECREATE");
    TDirectory* summary_dir = output_root.mkdir("summary_stats");
    if (summary_dir) {
        for (const auto& [scan_layer, layer_series] : layer_series_by_scan) {
            summary_dir->cd();
            std::string layer_dir_name = "layer" + std::to_string(scan_layer);
            TDirectory* layer_dir = summary_dir->mkdir(layer_dir_name.c_str());
            if (!layer_dir) {
                continue;
            }
            TDirectory* efficiencies_dir = layer_dir->mkdir("efficiencies");
            TDirectory* clusters_dir = layer_dir->mkdir("clusters");
            TDirectory* noise_dir = layer_dir->mkdir("noise");

            auto pick_directory = [&](std::string_view metric_name) -> TDirectory* {
                if (metric_name.rfind("eff_", 0) == 0 || metric_name.rfind("track_eff_", 0) == 0) {
                    return efficiencies_dir ? efficiencies_dir : layer_dir;
                }
                if (metric_name.rfind("avg_cluster_size", 0) == 0) {
                    return clusters_dir ? clusters_dir : layer_dir;
                }
                if (metric_name.rfind("noise_rate", 0) == 0) {
                    return noise_dir ? noise_dir : layer_dir;
                }
                return layer_dir;
            };

            for (const auto& [metric_name, series] : layer_series) {
                if (TDirectory* target_dir = pick_directory(metric_name)) {
                    target_dir->cd();
                }
                TMultiGraph* graph = new TMultiGraph();
                graph->SetName(metric_name.c_str());
                graph->SetTitle((metric_name + ";HV;Value").c_str());

                for (int layer = 0; layer < 3; ++layer) {
                    const auto& xs = series.x[layer];
                    const auto& ys = series.y[layer];
                    if (xs.empty()) { continue; }
                    TGraph* layer_graph = new TGraph(static_cast<int>(xs.size()), xs.data(), ys.data());
                    std::string layer_name = metric_name + "_layer" + std::to_string(layer);
                    layer_graph->SetName(layer_name.c_str());
                    layer_graph->SetMarkerStyle(20 + layer);
                    layer_graph->SetMarkerColor(1 + layer);
                    layer_graph->SetLineColor(1 + layer);
                    graph->Add(layer_graph, "P");
                }

                graph->Write();
            }

            const auto& scalar_x = scalar_x_by_scan[scan_layer];
            const auto& scalar_y = scalar_y_by_scan[scan_layer];
            for (const auto& [metric_name, xs] : scalar_x) {
                const auto& ys = scalar_y.at(metric_name);
                if (xs.empty()) { continue; }
                if (TDirectory* target_dir = pick_directory(metric_name)) {
                    target_dir->cd();
                }
                TGraph* graph = new TGraph(static_cast<int>(xs.size()), xs.data(), ys.data());
                graph->SetName(metric_name.c_str());
                graph->SetTitle((metric_name + ";HV;Value").c_str());
                graph->SetMarkerStyle(20);
                graph->Write();
            }
        }
    }

    output_root.Write();
    output_root.Close();

    summary_root->Close();
    delete summary_root;
}

// Accessors
const TGraph* DataPlotter::getGraphForMetric(const std::string& metric_name) const {
    std::filesystem::path config_stem = std::filesystem::path(config_path).stem();
    std::filesystem::path analysis_root_path = output_directory / (config_stem.string() + "_analysis.root");
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

// Utility functions
void DataPlotter::removeLayerFromGraph(TGraph* graph, int layer_to_remove) {
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