#include <iostream>
#include <algorithm>
#include <unordered_set>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "configParser.hpp"
#include "dataPlotter.hpp"
#include "plotStyler.hpp"
#include "plotBatchExporter.hpp"

#include "utils.hpp"
#include "constants.hpp"

#include "core/hit.hpp"
#include "core/event.hpp"
#include "core/cluster.hpp"
#include "core/track.hpp"
#include "analysis/dataAnalyzer.hpp"

// ==========================================================================================
// Analysis utility/helper namespaces for plotting and calculating statistics
// ==========================================================================================
namespace perFileHelpers {

// Anonymous namespace for internal helper functions for per-file plotting
namespace {
    int remapStrip(int rawStrip) {
        for (const auto& col : columnShifts)
            if (rawStrip >= col.start && rawStrip <= col.end)
                return rawStrip + col.shift;
        return rawStrip;
    }
}

void plotStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* strip_dir = analysis_dir->GetDirectory("strip");
    strip_dir->cd();

    TTree* input_data_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    if (!input_data_tree || input_data_tree->IsZombie()) {
        std::cerr << "Error: Invalid input data tree for analysis." << std::endl;
        return;
    }
    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create strip distributions for each layer and side
    TTreeReader readerInputData(input_data_tree);
    TTreeReader readerProcData(proc_tree);
    TTreeReaderValue<std::vector<int>> raw_time1(readerInputData, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>> raw_time2(readerInputData, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");

    // Create histograms using arrays
    const int nConfigs = 2;
    const char* categories[nConfigs] = {
        "strip_eta1", "strip_eta2"
    };

    std::map<std::string, std::map<int, TH1*>> strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH1F(Form("h1d_%s_layer%d", categories[c], layer),
                            Form("Layer %d;Strip;Hits", layer),
                            24, 0, 24);
            strip_histograms[categories[c]][layer] = hist;
        }
    }

    while (readerInputData.Next() && readerProcData.Next()) {
        for (size_t i = 0; i < strips->size(); ++i) {
            int layer = (*layers)[i];
            int strip = remapStrip((*strips)[i]);
            
            if ((*raw_time1)[i] != 0) strip_histograms["strip_eta1"][layer]->Fill(strip);
            if ((*raw_time2)[i] != 0) strip_histograms["strip_eta2"][layer]->Fill(strip);
        }
    }

    // Write histograms to file and clean up
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            strip_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete strip_histograms[categories[c]][layer];
        }
    }
}

void plotDtVsStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* dt_strip_dir = analysis_dir->GetDirectory("dt_strip");
    dt_strip_dir->cd();

    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }
    if (!track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid track reconstruction tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create 2D histograms for dt vs strip for each layer and valid track category
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> dts(readerProcData, "proc_dt_time1_time2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    // Create histograms using arrays
    const int nConfigs = 3;
    const char* suffixes[nConfigs] = {"dt_strip_all", "dt_strip_valid1", "dt_strip_valid2"};

    std::map<std::string, std::map<int, TH2*>> dt_strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(Form("h2d_%s_layer%d", suffixes[c], layer),
                            Form("Layer %d;Strip;dt [Ticks]; Entries", layer),
                            24, 0, 24, 40, -20, 20);
            dt_strip_histograms[suffixes[c]][layer] = hist;
        }
    }

    while (readerProcData.Next() && readerTrackData.Next()) {
        for (size_t i = 0; i < strips->size(); ++i) {
            int layer = (*layers)[i];
            int strip = remapStrip((*strips)[i]);
            int dt = (*dts)[i];
            
            dt_strip_histograms["dt_strip_all"][layer]->Fill(strip, dt);
            if ((*in_valid_track_eta1)[i]) {
                dt_strip_histograms["dt_strip_valid1"][layer]->Fill(strip, dt);
            }
            if ((*in_valid_track_eta2)[i]) {
                dt_strip_histograms["dt_strip_valid2"][layer]->Fill(strip, dt);
            }
        }
    }

    // Write histograms to file and clean up
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            dt_strip_histograms[suffixes[c]][layer]->Write("", TObject::kOverwrite);
            delete dt_strip_histograms[suffixes[c]][layer];
        }
    }
}

void plotToTVsStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* tot_strip_dir = analysis_dir->GetDirectory("tot_strip");
    tot_strip_dir->cd();

    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }
    if (!track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid track reconstruction tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create 2D histograms for dt vs strip for each layer and valid track category
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> proc_tot1(readerProcData, "proc_tot1");
    TTreeReaderValue<std::vector<int>> proc_tot2(readerProcData, "proc_tot2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    // Create histograms using arrays
    const int nConfigs = 6;
    const char* suffixes[nConfigs] = {"tot1_strip_all", "tot2_strip_all",
        "tot1_strip_valid1", "tot2_strip_valid1", "tot1_strip_valid2", "tot2_strip_valid2"};

    std::map<std::string, std::map<int, TH2*>> tot_strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(Form("h2d_%s_layer%d", suffixes[c], layer),
                            Form("Layer %d;Strip;ToT [Ticks]; Entries", layer),
                            24, 0, 24, 100, 0, 100);
            tot_strip_histograms[suffixes[c]][layer] = hist;
        }
    }

    while (readerProcData.Next() && readerTrackData.Next()) {
        for (size_t i = 0; i < strips->size(); ++i) {
            int layer = (*layers)[i];
            
            int strip = remapStrip((*strips)[i]);
            int tot1 = (*proc_tot1)[i];
            int tot2 = (*proc_tot2)[i];

            
            tot_strip_histograms["tot1_strip_all"][layer]->Fill(strip, tot1);
            tot_strip_histograms["tot2_strip_all"][layer]->Fill(strip, tot2);
            if ((*in_valid_track_eta1)[i]) {
                tot_strip_histograms["tot1_strip_valid1"][layer]->Fill(strip, tot1);
                tot_strip_histograms["tot2_strip_valid1"][layer]->Fill(strip, tot2);
            } else if ((*in_valid_track_eta2)[i]) {
                tot_strip_histograms["tot1_strip_valid2"][layer]->Fill(strip, tot1);
                tot_strip_histograms["tot2_strip_valid2"][layer]->Fill(strip, tot2);
            }
        }
    }

    // Write histograms to file and clean up
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {   
            tot_strip_histograms[suffixes[c]][layer]->Write("", TObject::kOverwrite);
            delete tot_strip_histograms[suffixes[c]][layer];
        }
    }

}

void plotMultiplicityAndDelayVsStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* mult_strip_dir = analysis_dir->GetDirectory("multiplicity_strip");
    TDirectory* delay_strip_dir = analysis_dir->GetDirectory("delay_strip");

    TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }
    if (!track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid track reconstruction tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create 2D histograms for dt vs strip for each layer and valid track category
    TTreeReader readerInputData(input_tree);
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> raw_time1(readerInputData, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>> raw_time2(readerInputData, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>> rise(readerInputData, "hit_rise");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> time1(readerProcData, "proc_time1");
    TTreeReaderValue<std::vector<int>> time2(readerProcData, "proc_time2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    // Create histograms using arrays
    const int nConfigs = 6;
    const char* categories[nConfigs] = {
        "mult1_all", "mult2_all",           // Multiplicity for all hits based on time1 and time2
        "mult1_valid1", "mult2_valid1",     // Multiplicity for hits in valid tracks on eta1 side based on time1 and time2
        "mult1_valid2", "mult2_valid2"      // Multiplicity for hits in valid tracks on eta2 side based on time1 and time2
    };

    std::map<std::string, std::map<int, TH2*>> multiplicity_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(Form("h2d_%s_layer%d", categories[c], layer),
                            Form("Layer %d;Strip;Multiplicity; Entries", layer),
                            24, 0, 24, 9, 1, 10);
            multiplicity_histograms[categories[c]][layer] = hist;
        }
    }
    std::map<std::string, std::map<int, TH2*>> delay_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(Form("h2d_delay_%s_layer%d", categories[c], layer),
                            Form("Layer %d;Strip;Delay from First Hit [Ticks]; Entries", layer),
                            24, 0, 24, 100, 0, 100);
            delay_histograms[categories[c]][layer] = hist;
        }
    }

    while (readerInputData.Next() && readerProcData.Next() && readerTrackData.Next()) {

        // Create lookup tables for each entry
        std::map<int, std::map<int, int>> counts[6];    // layer -> strip -> count (for each of the 6 categories)
        std::map<int, std::map<int, std::vector<int>>> delays[6];

        for (size_t i = 0; i < strips->size(); ++i) {

            if ((*rise)[i] == 0) continue; // Skip falling edge hits

            int layer = (*layers)[i];
            int strip = remapStrip((*strips)[i]);

            if ((*raw_time1)[i] != 0) {
                counts[0][layer][strip]++;
                delays[0][layer][strip].push_back((*time1)[i]);
                if ((*in_valid_track_eta1)[i]) {
                    counts[2][layer][strip]++;
                    delays[2][layer][strip].push_back((*time1)[i]);
                }
                if ((*in_valid_track_eta2)[i]) {
                    counts[4][layer][strip]++;
                    delays[4][layer][strip].push_back((*time1)[i]);
                }
            }
            if ((*raw_time2)[i] != 0) {
                counts[1][layer][strip]++;
                delays[1][layer][strip].push_back((*time2)[i]);
                if ((*in_valid_track_eta1)[i]) {
                    counts[3][layer][strip]++;
                    delays[3][layer][strip].push_back((*time2)[i]);
                }
                if ((*in_valid_track_eta2)[i]) {
                    counts[5][layer][strip]++;
                    delays[5][layer][strip].push_back((*time2)[i]);
                }
            }
        }

        // Fill histograms for this entry based on the counts and delays in the lookup tables
        for (int c = 0; c < nConfigs; ++c) {
            for (int layer : {0, 1, 2}) {
                for (const auto& [strip, count] : counts[c][layer]) {
                    multiplicity_histograms[categories[c]][layer]->Fill(strip, count);
                }

                for (const auto& [strip, delay_vec] : delays[c][layer]) {
                    if (!delay_vec.empty()) {
                        auto min_delay = std::min_element(delay_vec.begin(), delay_vec.end());
                        size_t min_index = std::distance(delay_vec.begin(), min_delay);
                        for (size_t i = 0; i < delay_vec.size(); ++i) {
                            if (i == min_index) continue; // Skip the first hit (minimum delay)
                            delay_histograms[categories[c]][layer]->Fill(strip, delay_vec[i] - *min_delay);
                        }
                    }
                }
            }
        }
    }

    // Write histograms to file and clean up
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            mult_strip_dir->cd();
            multiplicity_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete multiplicity_histograms[categories[c]][layer];
            delay_strip_dir->cd();
            delay_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete delay_histograms[categories[c]][layer];
        }
    }
}

}

// ==========================================================================================
// DataAnalyzer class implementation for analyzing processed DCT data and plotting results
// ==========================================================================================
DataAnalyzer::DataAnalyzer(const std::string& config_file_path, const std::string& output_directory_path)
    : _config_path({config_file_path})
    , _output_directory(output_directory_path) {
}

DataAnalyzer::~DataAnalyzer() {
}

void DataAnalyzer::producePerFileStats(TFile* input_file) {

    // Set up output directories for per-file statistics and plots
    std::vector<std::string> dir_names = {
        "analysis/strip",
        "analysis/dt_strip",
        "analysis/tot_strip",
        "analysis/multiplicity_strip",
        "analysis/delay_strip"
    };
    PathUtils::setupDirectories(input_file, dir_names);

    // Produce relevant plots for this file using helper functions
    perFileHelpers::plotStrip(input_file);
    perFileHelpers::plotDtVsStrip(input_file);
    perFileHelpers::plotToTVsStrip(input_file);
    perFileHelpers::plotMultiplicityAndDelayVsStrip(input_file);

    // Write analysis directory into the file
    input_file->Write("", TObject::kOverwrite);

    // Export per-file plots to PDF
    std::string root_file_path = input_file->GetName();
    std::string target_plots_dir = (_output_directory / "plots" / std::filesystem::path(input_file->GetName()).stem()).string();
    PlotterHelpers::BatchExporter::autoExportToATLASPDF(root_file_path, target_plots_dir);
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

    // Metadata variables to be stored in the summary tree
    std::string entry_name;
    std::string measurement_type;

    std::string mixture;
    bool source;
    double filter;
    int lv_setting;

    int scanned_layer;
    double scanned_hv = 0.0;
    double other_hv = 0.0;

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

    double noise_rate_error = 0.0;
    double noise_rate_eta1_error[3] = {0.0, 0.0, 0.0};
    double noise_rate_eta2_error[3] = {0.0, 0.0, 0.0};

    // Both scans relevant statistics
    double avg_cluster_size_eta1 = 0.0;
    double avg_cluster_size_eta2 = 0.0;
    double avg_cluster_size_eta1_layers[3] = {0.0, 0.0, 0.0};
    double avg_cluster_size_eta2_layers[3] = {0.0, 0.0, 0.0};

    double avg_cluster_size_eta1_error = 0.0;
    double avg_cluster_size_eta2_error = 0.0;
    double avg_cluster_size_eta1_layers_error[3] = {0.0, 0.0, 0.0};
    double avg_cluster_size_eta2_layers_error[3] = {0.0, 0.0, 0.0};

    // Set up branches for the summary tree
    summary_tree->Branch("name", &entry_name);
    summary_tree->Branch("measurement_type", &measurement_type);
    summary_tree->Branch("mixture", &mixture);
    summary_tree->Branch("source", &source);
    summary_tree->Branch("filter", &filter);
    summary_tree->Branch("lv_setting", &lv_setting);
    summary_tree->Branch("scanned_layer", &scanned_layer);
    summary_tree->Branch("scanned_hv", &scanned_hv);
    summary_tree->Branch("other_hv", &other_hv);

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

    summary_tree->Branch("avg_cluster_size_eta1", &avg_cluster_size_eta1);
    summary_tree->Branch("avg_cluster_size_eta2", &avg_cluster_size_eta2);
    summary_tree->Branch("avg_cluster_size_eta1_layers", avg_cluster_size_eta1_layers, "avg_cluster_size_eta1_layers[3]/D");
    summary_tree->Branch("avg_cluster_size_eta2_layers", avg_cluster_size_eta2_layers, "avg_cluster_size_eta2_layers[3]/D");
    summary_tree->Branch("avg_cluster_size_eta1_error", &avg_cluster_size_eta1_error);
    summary_tree->Branch("avg_cluster_size_eta2_error", &avg_cluster_size_eta2_error);
    summary_tree->Branch("avg_cluster_size_eta1_layers_error", avg_cluster_size_eta1_layers_error, "avg_cluster_size_eta1_layers_error[3]/D");
    summary_tree->Branch("avg_cluster_size_eta2_layers_error", avg_cluster_size_eta2_layers_error, "avg_cluster_size_eta2_layers_error[3]/D");

    summary_tree->Branch("noise_rate", &noise_rate);
    summary_tree->Branch("noise_rate_eta1", noise_rate_eta1, "noise_rate_eta1[3]/D");
    summary_tree->Branch("noise_rate_eta2", noise_rate_eta2, "noise_rate_eta2[3]/D");
    summary_tree->Branch("noise_rate_error", &noise_rate_error, "noise_rate_error[3]/D");
    summary_tree->Branch("noise_rate_eta1_error", noise_rate_eta1_error, "noise_rate_eta1_error[3]/D");
    summary_tree->Branch("noise_rate_eta2_error", noise_rate_eta2_error, "noise_rate_eta2_error[3]/D");

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

        // Calculate and fill per-file relevant statistics for this measurement entry and save into the input ROOT file
        producePerFileStats(input_file);

        // Extract metadata values for this measurement entry
        entry_name = metadata_entry.name;
        measurement_type = metadata_entry.measurement_type;
        mixture = metadata_entry.mixture;
        source = metadata_entry.source;
        filter = metadata_entry.filter;
        lv_setting = metadata_entry.lv_setting;
        scanned_layer = metadata_entry.scanned_layer;
        scanned_hv = metadata_entry.scanned_hv;
        other_hv = metadata_entry.other_hv;

        // Extract relevant data/statistics from the input ROOT file for this measurement entry
        MeasurementData stats;

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

        // Calculate noise rates per layer per side using raw hit times from InputData
        TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
        TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
        if (input_tree && proc_tree) {
            TTreeReader readerInputData(input_tree);
            TTreeReader readerProcData(proc_tree);
            TTreeReaderValue<std::vector<int>> channels(readerInputData, "hit_channel");
            TTreeReaderValue<std::vector<int>> time1(readerProcData, "proc_time1");
            TTreeReaderValue<std::vector<int>> time2(readerProcData, "proc_time2");

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

            while (readerInputData.Next() && readerProcData.Next()) {
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
        scan.data.push_back(stats);

        input_file->Close();
        delete input_file;
    }

    summary_root_file.cd();
    summary_tree->Write("", TObject::kOverwrite);
    summary_root_file.Close();
}
