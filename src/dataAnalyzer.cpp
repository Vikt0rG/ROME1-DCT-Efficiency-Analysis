#include <algorithm>
#include <array>
#include <unordered_set>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "dataAnalyzer.hpp"
#include "configParser.hpp"

/// TODO:
// - Check ToT calculation
// - Update root files with file-specific stats (see below)
// - Add noise rate calculations for each separate strip
// - Add time resolation calculations for each layer and side 


// ==========================================================================================
// DataAnalyzer class implementation for analyzing processed DCT data and plotting results
// ==========================================================================================
DataAnalyzer::DataAnalyzer(const std::string& config_file_path)
    : config_path({config_file_path}) {
}

DataAnalyzer::~DataAnalyzer() {
}

void DataAnalyzer::producePerFileStats() {
    // Function for producing per-file relevant statistics for each measurement entry such as:
    // proc_tot vs proc_strip for all three layers, both sides
    // proc_dt vs proc_strip for all three layers, both sides
    // heatmap of proc_strip vs proc_dt for all three layers, both sides
    // heatmap of hit_multiplicity vs proc_strip for all three layers, both sides
    // heatmap of proc_dt of consecutive hits vs proc_strip for all three layers, both sides
    // time resolution distributions for each layer and side

}

void DataAnalyzer::produceSummaryStats() {

    // Clear any existing entries and summaries before processing new config files
    parsed_entries.clear();
    summaries.clear();

    // Process each config file and build the list of measurement entries and summaries
    std::filesystem::create_directories(output_directory);

    // Prepare a summary object/root file for this config
    SummaryStats summary;
    summary.config_path = config_path;
    summary.entries = ConfigUtils::parseMeasurementEntries(config_path);
    if (!summary.entries.empty()) {
        summary.measurement_type = summary.entries.front().measurement_type;
    }

    parsed_entries.insert(parsed_entries.end(), summary.entries.begin(), summary.entries.end());
    std::filesystem::path config_stem = std::filesystem::path(config_path).stem();
    std::filesystem::path summary_root_path = output_directory / (config_stem.string() + "_summary.root");

    // Create a ROOT file and tree to store summary statistics for this config
    TFile summary_root(summary_root_path.string().c_str(), "RECREATE");
    TTree* summary_tree = new TTree("summary", "summary");

    // Metadata variables to be stored in the summary tree
    std::string entry_name;
    int layer;
    double hv = 0.0;

    // Efficiency scan relevant statistics
    double eta1_efficiency_external_summary[3] = {0.0, 0.0, 0.0};
    double eta2_efficiency_external_summary[3] = {0.0, 0.0, 0.0};
    double eta_or_efficiency_external_summary[3] = {0.0, 0.0, 0.0};
    double eta_and_efficiency_external_summary[3] = {0.0, 0.0, 0.0};

    double eta1_efficiency_external_error_summary[3] = {0.0, 0.0, 0.0};
    double eta2_efficiency_external_error_summary[3] = {0.0, 0.0, 0.0};
    double eta_or_efficiency_external_error_summary[3] = {0.0, 0.0, 0.0};
    double eta_and_efficiency_external_error_summary[3] = {0.0, 0.0, 0.0};

    double eta1_efficiency_rpc_summary[3] = {0.0, 0.0, 0.0};
    double eta2_efficiency_rpc_summary[3] = {0.0, 0.0, 0.0};
    double eta_or_efficiency_rpc_summary[3] = {0.0, 0.0, 0.0};
    double eta_and_efficiency_rpc_summary[3] = {0.0, 0.0, 0.0};

    double eta1_efficiency_rpc_error_summary[3] = {0.0, 0.0, 0.0};
    double eta2_efficiency_rpc_error_summary[3] = {0.0, 0.0, 0.0};
    double eta_or_efficiency_rpc_error_summary[3] = {0.0, 0.0, 0.0};
    double eta_and_efficiency_rpc_error_summary[3] = {0.0, 0.0, 0.0};

    double track_eta1_efficiency_external_summary[3] = {0.0, 0.0, 0.0};
    double track_eta2_efficiency_external_summary[3] = {0.0, 0.0, 0.0};
    double track_eta_or_efficiency_external_summary[3] = {0.0, 0.0, 0.0};
    double track_eta_and_efficiency_external_summary[3] = {0.0, 0.0, 0.0};

    double track_eta1_efficiency_external_error_summary[3] = {0.0, 0.0, 0.0};
    double track_eta2_efficiency_external_error_summary[3] = {0.0, 0.0, 0.0};
    double track_eta_or_efficiency_external_error_summary[3] = {0.0, 0.0, 0.0};
    double track_eta_and_efficiency_external_error_summary[3] = {0.0, 0.0, 0.0};

    double track_eta1_efficiency_rpc_summary[3] = {0.0, 0.0, 0.0};
    double track_eta2_efficiency_rpc_summary[3] = {0.0, 0.0, 0.0};
    double track_eta_or_efficiency_rpc_summary[3] = {0.0, 0.0, 0.0};
    double track_eta_and_efficiency_rpc_summary[3] = {0.0, 0.0, 0.0};

    double track_eta1_efficiency_rpc_error_summary[3] = {0.0, 0.0, 0.0};
    double track_eta2_efficiency_rpc_error_summary[3] = {0.0, 0.0, 0.0};
    double track_eta_or_efficiency_rpc_error_summary[3] = {0.0, 0.0, 0.0};
    double track_eta_and_efficiency_rpc_error_summary[3] = {0.0, 0.0, 0.0};

    // Noise scan relevant statistics
    double noise_rate = 0.0;
    double noise_rate_eta1[3] = {0.0, 0.0, 0.0};
    double noise_rate_eta2[3] = {0.0, 0.0, 0.0};

    // Both scans relevant statistics
    double avg_cluster_size_eta2 = 0.0;
    double avg_cluster_size_eta1 = 0.0;
    double avg_cluster_size_eta1_layers[3] = {0.0, 0.0, 0.0};
    double avg_cluster_size_eta2_layers[3] = {0.0, 0.0, 0.0};

    // Set up branches for the summary tree
    summary_tree->Branch("name", &entry_name);
    summary_tree->Branch("layer", &layer);
    summary_tree->Branch("hv", &hv);
    summary_tree->Branch("avg_cluster_size_eta1", &avg_cluster_size_eta1);
    summary_tree->Branch("avg_cluster_size_eta2", &avg_cluster_size_eta2);
    summary_tree->Branch("avg_cluster_size_eta1_layers", avg_cluster_size_eta1_layers, "avg_cluster_size_eta1_layers[3]/D");
    summary_tree->Branch("avg_cluster_size_eta2_layers", avg_cluster_size_eta2_layers, "avg_cluster_size_eta2_layers[3]/D");

    // Save efficiencies as graphs??
    summary_tree->Branch("eff_eta1_external", eta1_efficiency_external_summary, "eff_eta1_external[3]/D");
    summary_tree->Branch("eff_eta2_external", eta2_efficiency_external_summary, "eff_eta2_external[3]/D");
    summary_tree->Branch("eff_or_external", eta_or_efficiency_external_summary, "eff_or_external[3]/D");
    summary_tree->Branch("eff_and_external", eta_and_efficiency_external_summary, "eff_and_external[3]/D");
    summary_tree->Branch("eff_eta1_external_error", eta1_efficiency_external_error_summary, "eff_eta1_external_error[3]/D");
    summary_tree->Branch("eff_eta2_external_error", eta2_efficiency_external_error_summary, "eff_eta2_external_error[3]/D");
    summary_tree->Branch("eff_or_external_error", eta_or_efficiency_external_error_summary, "eff_or_external_error[3]/D");
    summary_tree->Branch("eff_and_external_error", eta_and_efficiency_external_error_summary, "eff_and_external_error[3]/D");

    summary_tree->Branch("eff_eta1_rpc", eta1_efficiency_rpc_summary, "eff_eta1_rpc[3]/D");
    summary_tree->Branch("eff_eta2_rpc", eta2_efficiency_rpc_summary, "eff_eta2_rpc[3]/D");
    summary_tree->Branch("eff_or_rpc", eta_or_efficiency_rpc_summary, "eff_or_rpc[3]/D");
    summary_tree->Branch("eff_and_rpc", eta_and_efficiency_rpc_summary, "eff_and_rpc[3]/D");
    summary_tree->Branch("eff_eta1_rpc_error", eta1_efficiency_rpc_error_summary, "eff_eta1_rpc_error[3]/D");
    summary_tree->Branch("eff_eta2_rpc_error", eta2_efficiency_rpc_error_summary, "eff_eta2_rpc_error[3]/D");
    summary_tree->Branch("eff_or_rpc_error", eta_or_efficiency_rpc_error_summary, "eff_or_rpc_error[3]/D");
    summary_tree->Branch("eff_and_rpc_error", eta_and_efficiency_rpc_error_summary, "eff_and_rpc_error[3]/D");

    summary_tree->Branch("track_eff_eta1_external", track_eta1_efficiency_external_summary, "track_eff_eta1_external[3]/D");
    summary_tree->Branch("track_eff_eta2_external", track_eta2_efficiency_external_summary, "track_eff_eta2_external[3]/D");
    summary_tree->Branch("track_eff_or_external", track_eta_or_efficiency_external_summary, "track_eff_or_external[3]/D");
    summary_tree->Branch("track_eff_and_external", track_eta_and_efficiency_external_summary, "track_eff_and_external[3]/D");
    summary_tree->Branch("track_eff_eta1_external_error", track_eta1_efficiency_external_error_summary, "track_eff_eta1_external_error[3]/D");
    summary_tree->Branch("track_eff_eta2_external_error", track_eta2_efficiency_external_error_summary, "track_eff_eta2_external_error[3]/D");
    summary_tree->Branch("track_eff_or_external_error", track_eta_or_efficiency_external_error_summary, "track_eff_or_external_error[3]/D");
    summary_tree->Branch("track_eff_and_external_error", track_eta_and_efficiency_external_error_summary, "track_eff_and_external_error[3]/D");

    summary_tree->Branch("track_eff_eta1_rpc", track_eta1_efficiency_rpc_summary, "track_eff_eta1_rpc[3]/D");
    summary_tree->Branch("track_eff_eta2_rpc", track_eta2_efficiency_rpc_summary, "track_eff_eta2_rpc[3]/D");
    summary_tree->Branch("track_eff_or_rpc", track_eta_or_efficiency_rpc_summary, "track_eff_or_rpc[3]/D");
    summary_tree->Branch("track_eff_and_rpc", track_eta_and_efficiency_rpc_summary, "track_eff_and_rpc[3]/D");
    summary_tree->Branch("track_eff_eta1_rpc_error", track_eta1_efficiency_rpc_error_summary, "track_eff_eta1_rpc_error[3]/D");
    summary_tree->Branch("track_eff_eta2_rpc_error", track_eta2_efficiency_rpc_error_summary, "track_eff_eta2_rpc_error[3]/D");
    summary_tree->Branch("track_eff_or_rpc_error", track_eta_or_efficiency_rpc_error_summary, "track_eff_or_rpc_error[3]/D");
    summary_tree->Branch("track_eff_and_rpc_error", track_eta_and_efficiency_rpc_error_summary, "track_eff_and_rpc_error[3]/D");

    summary_tree->Branch("noise_rate", &noise_rate);
    summary_tree->Branch("noise_rate_eta1", noise_rate_eta1, "noise_rate_eta1[3]/D");
    summary_tree->Branch("noise_rate_eta2", noise_rate_eta2, "noise_rate_eta2[3]/D");

    // Process each measurement entry in this config file
    for (const auto& entry : summary.entries) {
        if (entry.root_file.empty()) {
            continue;
        }

        TFile* input_file = TFile::Open(entry.root_file.c_str(), "READ");
        if (!input_file || input_file->IsZombie()) {
            if (input_file) {
                input_file->Close();
                delete input_file;
            }
            continue;
        }

        // TODO: makeAnalysisTree(&input_file)

        // Extract relevant statistics from the input ROOT file for this measurement entry
        PerFileStats stats;
        stats.name = entry.name;
        stats.layer = entry.layer;
        stats.hv = entry.hv;

        // Extract efficiency values from histograms for this measurement entry
        for (int layer = 0; layer < 3; ++layer) {

            // External trigger only
            std::string hist_path = "efficiencies_histograms/external_trigger/eff_external_trigger_layer" + std::to_string(layer);
            TH1* hist = dynamic_cast<TH1*>(input_file->Get(hist_path.c_str()));
            if (!hist) {
                continue;
            }
            stats.efficiency_results.eta1_efficiency_external[layer] = hist->GetBinContent(1);
            stats.efficiency_results.eta2_efficiency_external[layer] = hist->GetBinContent(2);
            stats.efficiency_results.eta_or_efficiency_external[layer] = hist->GetBinContent(3);
            stats.efficiency_results.eta_and_efficiency_external[layer] = hist->GetBinContent(4);

            stats.efficiency_results.eta1_efficiency_external_error[layer] = hist->GetBinError(1);
            stats.efficiency_results.eta2_efficiency_external_error[layer] = hist->GetBinError(2);
            stats.efficiency_results.eta_or_efficiency_external_error[layer] = hist->GetBinError(3);
            stats.efficiency_results.eta_and_efficiency_external_error[layer] = hist->GetBinError(4);

            // External trigger + RPC
            std::string hist_path_rpc = "efficiencies_histograms/external_plus_rpc_trigger/eff_rpc_layer" + std::to_string(layer);
            TH1* hist_rpc = dynamic_cast<TH1*>(input_file->Get(hist_path_rpc.c_str()));
            if (!hist_rpc) {
                continue;
            }
            stats.efficiency_results.eta1_efficiency_rpc[layer] = hist_rpc->GetBinContent(1);
            stats.efficiency_results.eta2_efficiency_rpc[layer] = hist_rpc->GetBinContent(2);
            stats.efficiency_results.eta_or_efficiency_rpc[layer] = hist_rpc->GetBinContent(3);
            stats.efficiency_results.eta_and_efficiency_rpc[layer] = hist_rpc->GetBinContent(4);

            stats.efficiency_results.eta1_efficiency_rpc_error[layer] = hist_rpc->GetBinError(1);
            stats.efficiency_results.eta2_efficiency_rpc_error[layer] = hist_rpc->GetBinError(2);
            stats.efficiency_results.eta_or_efficiency_rpc_error[layer] = hist_rpc->GetBinError(3);
            stats.efficiency_results.eta_and_efficiency_rpc_error[layer] = hist_rpc->GetBinError(4);

            // Track external trigger only
            std::string hist_path_track = "efficiencies_histograms/track_external_trigger/track_eff_external_trigger_layer" + std::to_string(layer);
            TH1* hist_track = dynamic_cast<TH1*>(input_file->Get(hist_path_track.c_str()));
            if (!hist_track) {
                continue;
            }
            stats.efficiency_results_tracks.track_eta1_efficiency_external[layer] = hist_track->GetBinContent(1);
            stats.efficiency_results_tracks.track_eta2_efficiency_external[layer] = hist_track->GetBinContent(2);
            stats.efficiency_results_tracks.track_eta_or_efficiency_external[layer] = hist_track->GetBinContent(3);
            stats.efficiency_results_tracks.track_eta_and_efficiency_external[layer] = hist_track->GetBinContent(4);

            stats.efficiency_results_tracks.track_eta1_efficiency_external_error[layer] = hist_track->GetBinError(1);
            stats.efficiency_results_tracks.track_eta2_efficiency_external_error[layer] = hist_track->GetBinError(2);
            stats.efficiency_results_tracks.track_eta_or_efficiency_external_error[layer] = hist_track->GetBinError(3);
            stats.efficiency_results_tracks.track_eta_and_efficiency_external_error[layer] = hist_track->GetBinError(4);

            // Track external trigger + RPC
            std::string hist_path_track_rpc = "efficiencies_histograms/track_external_plus_rpc_trigger/track_eff_rpc_layer" + std::to_string(layer);
            TH1* hist_track_rpc = dynamic_cast<TH1*>(input_file->Get(hist_path_track_rpc.c_str()));
            if (!hist_track_rpc) {
                continue;
            }
            stats.efficiency_results_tracks.track_eta1_efficiency_rpc[layer] = hist_track_rpc->GetBinContent(1);
            stats.efficiency_results_tracks.track_eta2_efficiency_rpc[layer] = hist_track_rpc->GetBinContent(2);
            stats.efficiency_results_tracks.track_eta_or_efficiency_rpc[layer] = hist_track_rpc->GetBinContent(3);
            stats.efficiency_results_tracks.track_eta_and_efficiency_rpc[layer] = hist_track_rpc->GetBinContent(4);

            stats.efficiency_results_tracks.track_eta1_efficiency_rpc_error[layer] = hist_track_rpc->GetBinError(1);
            stats.efficiency_results_tracks.track_eta2_efficiency_rpc_error[layer] = hist_track_rpc->GetBinError(2);
            stats.efficiency_results_tracks.track_eta_or_efficiency_rpc_error[layer] = hist_track_rpc->GetBinError(3);
            stats.efficiency_results_tracks.track_eta_and_efficiency_rpc_error[layer] = hist_track_rpc->GetBinError(4);
        }


        // Extract clusterization statistics from the "Clusterization" tree for this measurement entry
        TTree* cluster_tree = dynamic_cast<TTree*>(input_file->Get("Clusterization"));
        if (cluster_tree) {
            TTreeReader reader(cluster_tree);
            TTreeReaderValue<std::vector<int>> cluster_eta1(reader, "cluster_size_eta1");
            TTreeReaderValue<std::vector<int>> cluster_eta2(reader, "cluster_size_eta2");
            TTreeReaderValue<std::vector<int>> cluster_eta1_layer0(reader, "cluster_size_eta1_layer0");
            TTreeReaderValue<std::vector<int>> cluster_eta1_layer1(reader, "cluster_size_eta1_layer1");
            TTreeReaderValue<std::vector<int>> cluster_eta1_layer2(reader, "cluster_size_eta1_layer2");
            TTreeReaderValue<std::vector<int>> cluster_eta2_layer0(reader, "cluster_size_eta2_layer0");
            TTreeReaderValue<std::vector<int>> cluster_eta2_layer1(reader, "cluster_size_eta2_layer1");
            TTreeReaderValue<std::vector<int>> cluster_eta2_layer2(reader, "cluster_size_eta2_layer2");

            long long total_clusters_eta1 = 0;
            long long total_clusters_eta2 = 0;
            long long total_cluster_size_eta1 = 0;
            long long total_cluster_size_eta2 = 0;
            long long total_clusters_eta1_layers[3] = {0, 0, 0};
            long long total_clusters_eta2_layers[3] = {0, 0, 0};
            long long total_cluster_size_eta1_layers[3] = {0, 0, 0};
            long long total_cluster_size_eta2_layers[3] = {0, 0, 0};

            while (reader.Next()) {
                if (cluster_eta1.GetSetupStatus() == 0) {
                    total_clusters_eta1 += static_cast<long long>(cluster_eta1->size());
                    for (int size : *cluster_eta1) {
                        total_cluster_size_eta1 += size;
                    }
                }
                if (cluster_eta2.GetSetupStatus() == 0) {
                    total_clusters_eta2 += static_cast<long long>(cluster_eta2->size());
                    for (int size : *cluster_eta2) {
                        total_cluster_size_eta2 += size;
                    }
                }
                if (cluster_eta1_layer0.GetSetupStatus() == 0) {
                    total_clusters_eta1_layers[0] += static_cast<long long>(cluster_eta1_layer0->size());
                    for (int size : *cluster_eta1_layer0) {
                        total_cluster_size_eta1_layers[0] += size;
                    }
                }
                if (cluster_eta1_layer1.GetSetupStatus() == 0) {
                    total_clusters_eta1_layers[1] += static_cast<long long>(cluster_eta1_layer1->size());
                    for (int size : *cluster_eta1_layer1) {
                        total_cluster_size_eta1_layers[1] += size;
                    }
                }
                if (cluster_eta1_layer2.GetSetupStatus() == 0) {
                    total_clusters_eta1_layers[2] += static_cast<long long>(cluster_eta1_layer2->size());
                    for (int size : *cluster_eta1_layer2) {
                        total_cluster_size_eta1_layers[2] += size;
                    }
                }
                if (cluster_eta2_layer0.GetSetupStatus() == 0) {
                    total_clusters_eta2_layers[0] += static_cast<long long>(cluster_eta2_layer0->size());
                    for (int size : *cluster_eta2_layer0) {
                        total_cluster_size_eta2_layers[0] += size;
                    }
                }
                if (cluster_eta2_layer1.GetSetupStatus() == 0) {
                    total_clusters_eta2_layers[1] += static_cast<long long>(cluster_eta2_layer1->size());
                    for (int size : *cluster_eta2_layer1) {
                        total_cluster_size_eta2_layers[1] += size;
                    }
                }
                if (cluster_eta2_layer2.GetSetupStatus() == 0) {
                    total_clusters_eta2_layers[2] += static_cast<long long>(cluster_eta2_layer2->size());
                    for (int size : *cluster_eta2_layer2) {
                        total_cluster_size_eta2_layers[2] += size;
                    }
                }
            }

            if (total_clusters_eta1 > 0) {
                stats.avg_cluster_size_eta1 = static_cast<double>(total_cluster_size_eta1) / total_clusters_eta1;
            }
            if (total_clusters_eta2 > 0) {
                stats.avg_cluster_size_eta2 = static_cast<double>(total_cluster_size_eta2) / total_clusters_eta2;
            }
            for (int layer_idx = 0; layer_idx < 3; ++layer_idx) {
                if (total_clusters_eta1_layers[layer_idx] > 0) {
                    stats.avg_cluster_size_eta1_layers[layer_idx] =
                        static_cast<double>(total_cluster_size_eta1_layers[layer_idx]) / total_clusters_eta1_layers[layer_idx];
                }
                if (total_clusters_eta2_layers[layer_idx] > 0) {
                    stats.avg_cluster_size_eta2_layers[layer_idx] =
                        static_cast<double>(total_cluster_size_eta2_layers[layer_idx]) / total_clusters_eta2_layers[layer_idx];
                }
            }
        }

        // Calculate noise rate
        TTree* processed_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
        if (processed_tree) {
            TTreeReader reader(processed_tree);

            // Get number of total hits, number of events, and unique strips/layers in one pass.
            TTreeReaderValue<int> hits(reader, "n_hits");
            TTreeReaderValue<std::vector<int>> strips(reader, "proc_strip");
            TTreeReaderValue<std::vector<int>> layers(reader, "proc_layer");

            long long total_hits = 0;
            long long event_count = 0;
            std::unordered_set<int> unique_strips;
            std::unordered_set<int> unique_layers;

            while (reader.Next()) {
                total_hits += static_cast<long long>(*hits);
                event_count += 1;

                if (strips.GetSetupStatus() == 0 && layers.GetSetupStatus() == 0) {
                    for (int strip : *strips) {
                        unique_strips.insert(strip);
                    }
                    for (int layer : *layers) {
                        unique_layers.insert(layer);
                    }
                }
            }

            int n_strips = static_cast<int>(unique_strips.size());
            int n_layers = static_cast<int>(unique_layers.size());

            // Calculate the noise rate
            if (event_count > 0 && n_strips > 0) {
                stats.noise_rate = static_cast<double>(total_hits) / (event_count * TRIGGER_TIME_WINDOW * n_layers * n_strips * STRIP_WIDTH_CM * DETECTOR_LENGTH_CM);
            }
        }

        // Calculate split noise rates using raw hit times from InputData
        TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
        if (input_tree) {
            TTreeReader reader(input_tree);
            TTreeReaderValue<std::vector<int>> channels(reader, "hit_channel");
            TTreeReaderValue<std::vector<int>> time1(reader, "hit_time1");
            TTreeReaderValue<std::vector<int>> time2(reader, "hit_time2");

            std::array<long long, 3> total_hits_eta1 = {0, 0, 0};
            std::array<long long, 3> total_hits_eta2 = {0, 0, 0};
            std::array<std::unordered_set<int>, 3> unique_strips_eta1;
            std::array<std::unordered_set<int>, 3> unique_strips_eta2;

            auto channel_to_layer_strip = [](int channel, int& out_layer, int& out_strip) {
                const int layer = (channel % 24) / 8;
                const int column = channel / 24;
                const int strip = 8 * column + channel % 8;
                out_layer = layer;
                out_strip = strip;
            };

            while (reader.Next()) {
                if (channels.GetSetupStatus() != 0 || time1.GetSetupStatus() != 0 || time2.GetSetupStatus() != 0) {
                    continue;
                }
                const size_t n = std::min({channels->size(), time1->size(), time2->size()});
                for (size_t i = 0; i < n; ++i) {
                    int layer_idx = 0;
                    int strip = 0;
                    channel_to_layer_strip((*channels)[i], layer_idx, strip);
                    if (layer_idx < 0 || layer_idx >= 3) {
                        continue;
                    }
                    if ((*time1)[i] != 0) {
                        total_hits_eta1[layer_idx] += 1;
                        unique_strips_eta1[layer_idx].insert(strip);
                    }
                    if ((*time2)[i] != 0) {
                        total_hits_eta2[layer_idx] += 1;
                        unique_strips_eta2[layer_idx].insert(strip);
                    }
                }
            }

            const long long event_count = processed_tree ? processed_tree->GetEntries() : 0;
            for (int layer_idx = 0; layer_idx < 3; ++layer_idx) {
                const int n_strips_eta1 = static_cast<int>(unique_strips_eta1[layer_idx].size());
                const int n_strips_eta2 = static_cast<int>(unique_strips_eta2[layer_idx].size());
                if (event_count > 0 && n_strips_eta1 > 0) {
                    stats.noise_rate_eta1[layer_idx] = static_cast<double>(total_hits_eta1[layer_idx]) /
                                                     (event_count * TRIGGER_TIME_WINDOW * n_strips_eta1 * STRIP_WIDTH_CM * DETECTOR_LENGTH_CM);
                }
                if (event_count > 0 && n_strips_eta2 > 0) {
                    stats.noise_rate_eta2[layer_idx] = static_cast<double>(total_hits_eta2[layer_idx]) /
                                                     (event_count * TRIGGER_TIME_WINDOW * n_strips_eta2 * STRIP_WIDTH_CM * DETECTOR_LENGTH_CM);
                }
            }
        }

        entry_name = stats.name;
        layer = stats.layer;
        hv = stats.hv;
        avg_cluster_size_eta1 = stats.avg_cluster_size_eta1;
        avg_cluster_size_eta2 = stats.avg_cluster_size_eta2;
        for (int layer = 0; layer < 3; ++layer) {
            avg_cluster_size_eta1_layers[layer] = stats.avg_cluster_size_eta1_layers[layer];
            avg_cluster_size_eta2_layers[layer] = stats.avg_cluster_size_eta2_layers[layer];
            eta1_efficiency_external_summary[layer] = stats.efficiency_results.eta1_efficiency_external[layer];
            eta2_efficiency_external_summary[layer] = stats.efficiency_results.eta2_efficiency_external[layer];
            eta_or_efficiency_external_summary[layer] = stats.efficiency_results.eta_or_efficiency_external[layer];
            eta_and_efficiency_external_summary[layer] = stats.efficiency_results.eta_and_efficiency_external[layer];

            eta1_efficiency_external_error_summary[layer] = stats.efficiency_results.eta1_efficiency_external_error[layer];
            eta2_efficiency_external_error_summary[layer] = stats.efficiency_results.eta2_efficiency_external_error[layer];
            eta_or_efficiency_external_error_summary[layer] = stats.efficiency_results.eta_or_efficiency_external_error[layer];
            eta_and_efficiency_external_error_summary[layer] = stats.efficiency_results.eta_and_efficiency_external_error[layer];

            eta1_efficiency_rpc_summary[layer] = stats.efficiency_results.eta1_efficiency_rpc[layer];
            eta2_efficiency_rpc_summary[layer] = stats.efficiency_results.eta2_efficiency_rpc[layer];
            eta_or_efficiency_rpc_summary[layer] = stats.efficiency_results.eta_or_efficiency_rpc[layer];
            eta_and_efficiency_rpc_summary[layer] = stats.efficiency_results.eta_and_efficiency_rpc[layer];

            eta1_efficiency_rpc_error_summary[layer] = stats.efficiency_results.eta1_efficiency_rpc_error[layer];
            eta2_efficiency_rpc_error_summary[layer] = stats.efficiency_results.eta2_efficiency_rpc_error[layer];
            eta_or_efficiency_rpc_error_summary[layer] = stats.efficiency_results.eta_or_efficiency_rpc_error[layer];
            eta_and_efficiency_rpc_error_summary[layer] = stats.efficiency_results.eta_and_efficiency_rpc_error[layer];

            track_eta1_efficiency_external_summary[layer] = stats.efficiency_results_tracks.track_eta1_efficiency_external[layer];
            track_eta2_efficiency_external_summary[layer] = stats.efficiency_results_tracks.track_eta2_efficiency_external[layer];
            track_eta_or_efficiency_external_summary[layer] = stats.efficiency_results_tracks.track_eta_or_efficiency_external[layer];
            track_eta_and_efficiency_external_summary[layer] = stats.efficiency_results_tracks.track_eta_and_efficiency_external[layer];

            track_eta1_efficiency_external_error_summary[layer] = stats.efficiency_results_tracks.track_eta1_efficiency_external_error[layer];
            track_eta2_efficiency_external_error_summary[layer] = stats.efficiency_results_tracks.track_eta2_efficiency_external_error[layer];
            track_eta_or_efficiency_external_error_summary[layer] = stats.efficiency_results_tracks.track_eta_or_efficiency_external_error[layer];
            track_eta_and_efficiency_external_error_summary[layer] = stats.efficiency_results_tracks.track_eta_and_efficiency_external_error[layer];

            track_eta1_efficiency_rpc_summary[layer] = stats.efficiency_results_tracks.track_eta1_efficiency_rpc[layer];
            track_eta2_efficiency_rpc_summary[layer] = stats.efficiency_results_tracks.track_eta2_efficiency_rpc[layer];
            track_eta_or_efficiency_rpc_summary[layer] = stats.efficiency_results_tracks.track_eta_or_efficiency_rpc[layer];
            track_eta_and_efficiency_rpc_summary[layer] = stats.efficiency_results_tracks.track_eta_and_efficiency_rpc[layer];

            track_eta1_efficiency_rpc_error_summary[layer] = stats.efficiency_results_tracks.track_eta1_efficiency_rpc_error[layer];
            track_eta2_efficiency_rpc_error_summary[layer] = stats.efficiency_results_tracks.track_eta2_efficiency_rpc_error[layer];
            track_eta_or_efficiency_rpc_error_summary[layer] = stats.efficiency_results_tracks.track_eta_or_efficiency_rpc_error[layer];
            track_eta_and_efficiency_rpc_error_summary[layer] = stats.efficiency_results_tracks.track_eta_and_efficiency_rpc_error[layer];

            noise_rate_eta1[layer] = stats.noise_rate_eta1[layer];
            noise_rate_eta2[layer] = stats.noise_rate_eta2[layer];
        }
        noise_rate = stats.noise_rate;
        summary_tree->Fill();
        summary.per_file_stats.push_back(stats);

        input_file->Close();
        delete input_file;
    }

    summary_root.cd();
    summary_tree->Write();
    summary_root.Close();
    summaries.push_back(std::move(summary));

}