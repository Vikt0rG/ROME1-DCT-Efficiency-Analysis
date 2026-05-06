#include "dataAnalyzer.hpp"


// ==========================================================================================
// DataAnalyzer class implementation
// ==========================================================================================
/// Constructor and destructor for DataAnalyzer class
DataAnalyzer::DataAnalyzer() {
    output_file = nullptr;

    input_data_tree = nullptr;
    processed_data_tree = nullptr;
    clusterization_tree = nullptr;
    track_reconstruction_tree = nullptr;

    current_event = nullptr;
    current_event_number = 0;
    n_hits = 0;
    BC0 = 0;

    dt_max = -100;    // Default time window (ticks) - see constants.hpp
    dt_min = -180;    // Default time window (ticks) - see constants.hpp
    initializeCounters();
}

DataAnalyzer::~DataAnalyzer() {
    if (output_file) { 
        // Write all trees to file before closing
        if (input_data_tree) input_data_tree->Write();
        if (processed_data_tree) processed_data_tree->Write();
        if (clusterization_tree) clusterization_tree->Write();
        if (track_reconstruction_tree) track_reconstruction_tree->Write();
        output_file->Close(); 
        delete output_file; 
    }
}

/// Utility function to setup the output ROOT file and create the necessary trees
void DataAnalyzer::setupOutputFile() {
    output_file = new TFile("output.root", "RECREATE");
    input_data_tree = new TTree("InputData", "Raw hit data from the DCT");
    processed_data_tree = new TTree("ProcessedData", "Processed hit and event-level data");
    clusterization_tree = new TTree("Clusterization", "Cluster-level data");
    track_reconstruction_tree = new TTree("TrackReconstruction", "Track-level data");
}

/// TODO: Utility function to create input name from the input directory

/// Utility function to setup branches for all trees in the output file
void DataAnalyzer::setupBranches() {
    // Branch definitions for raw data tree
    input_data_tree->Branch("hit_clk", &hit_clk);
    input_data_tree->Branch("hit_channel", &hit_channel);
    input_data_tree->Branch("hit_raw_bcid", &hit_raw_bcid);
    input_data_tree->Branch("hit_bcid", &hit_bcid);
    input_data_tree->Branch("hit_time1", &hit_time1);
    input_data_tree->Branch("hit_time2", &hit_time2);
    input_data_tree->Branch("hit_rise", &hit_rise);
    input_data_tree->Branch("hit_raw_bcout", &hit_raw_bcout);
    input_data_tree->Branch("hit_bcout", &hit_bcout);

    // Branch definitions for processed data tree
    processed_data_tree->Branch("n_events", &current_event_number);
    processed_data_tree->Branch("n_hits", &n_hits);
    processed_data_tree->Branch("proc_layer", &proc_layer);
    processed_data_tree->Branch("proc_strip", &proc_strip);
    processed_data_tree->Branch("proc_time1", &proc_time1);
    processed_data_tree->Branch("proc_time2", &proc_time2);
    processed_data_tree->Branch("proc_dt_time1_time2", &proc_dt_time1_time2);
    processed_data_tree->Branch("proc_trigger_time", &proc_trigger_time);
    processed_data_tree->Branch("proc_dt_time1_trigger", &proc_dt_time1_trigger);
    processed_data_tree->Branch("proc_dt_time2_trigger", &proc_dt_time2_trigger);
    processed_data_tree->Branch("proc_tot1", &proc_tot1);
    processed_data_tree->Branch("proc_tot2", &proc_tot2);

    // Branch definitions for clusterization tree
    clusterization_tree->Branch("cluster_size_eta1", &cluster_size_eta1);
    clusterization_tree->Branch("cluster_size_eta2", &cluster_size_eta2);
    clusterization_tree->Branch("cluster_tot1", &cluster_tot1);
    clusterization_tree->Branch("cluster_tot2", &cluster_tot2);

    // Branch definitions for track reconstruction tree
    track_reconstruction_tree->Branch("track_length_eta1", &track_length_eta1);
    track_reconstruction_tree->Branch("track_length_eta2", &track_length_eta2);
    track_reconstruction_tree->Branch("track_width_eta1", &track_width_eta1);
    track_reconstruction_tree->Branch("track_width_eta2", &track_width_eta2);
    track_reconstruction_tree->Branch("track_size_eta1", &track_size_eta1);
    track_reconstruction_tree->Branch("track_size_eta2", &track_size_eta2);
}

/// Utility function to initialize efficiency counters and results
void DataAnalyzer::initializeCounters() {
    efficiency_counters = {};
    efficiency_counters_tracks = {};
    efficiency_results = {};
    efficiency_results_tracks = {};
}

/// Efficiency live calculation based on reconstructed tracks
void DataAnalyzer::updateEfficiencies() {
    // External trigger only efficiency results
    if (efficiency_counters.triggered_events_external > 0) {
        for (int layer = 0; layer < 3; layer++) {
            efficiency_results.eta1_efficiency_external[layer] = static_cast<double>(efficiency_counters.eta1_efficiency_counter[layer]) / efficiency_counters.triggered_events_external;
            efficiency_results.eta2_efficiency_external[layer] = static_cast<double>(efficiency_counters.eta2_efficiency_counter[layer]) / efficiency_counters.triggered_events_external;
            efficiency_results.eta_or_efficiency_external[layer] = static_cast<double>(efficiency_counters.eta_or_efficiency_counter[layer]) / efficiency_counters.triggered_events_external;
            efficiency_results.eta_and_efficiency_external[layer] = static_cast<double>(efficiency_counters.eta_and_efficiency_counter[layer]) / efficiency_counters.triggered_events_external;

            efficiency_results.eta1_efficiency_external_error[layer] = sqrt(efficiency_results.eta1_efficiency_external[layer] * (1 - efficiency_results.eta1_efficiency_external[layer]) / efficiency_counters.triggered_events_external);
            efficiency_results.eta2_efficiency_external_error[layer] = sqrt(efficiency_results.eta2_efficiency_external[layer] * (1 - efficiency_results.eta2_efficiency_external[layer]) / efficiency_counters.triggered_events_external);
            efficiency_results.eta_or_efficiency_external_error[layer] = sqrt(efficiency_results.eta_or_efficiency_external[layer] * (1 - efficiency_results.eta_or_efficiency_external[layer]) / efficiency_counters.triggered_events_external);
            efficiency_results.eta_and_efficiency_external_error[layer] = sqrt(efficiency_results.eta_and_efficiency_external[layer] * (1 - efficiency_results.eta_and_efficiency_external[layer]) / efficiency_counters.triggered_events_external);
        }
    }
    if (efficiency_counters_tracks.track_triggered_events_external > 0) {
        for (int layer = 0; layer < 3; layer++) {
            efficiency_results_tracks.track_eta1_efficiency_external[layer] = static_cast<double>(efficiency_counters_tracks.track_eta1_efficiency_counter[layer]) / efficiency_counters_tracks.track_triggered_events_external;
            efficiency_results_tracks.track_eta2_efficiency_external[layer] = static_cast<double>(efficiency_counters_tracks.track_eta2_efficiency_counter[layer]) / efficiency_counters_tracks.track_triggered_events_external;
            efficiency_results_tracks.track_eta_or_efficiency_external[layer] = static_cast<double>(efficiency_counters_tracks.track_eta_or_efficiency_counter[layer]) / efficiency_counters_tracks.track_triggered_events_external;
            efficiency_results_tracks.track_eta_and_efficiency_external[layer] = static_cast<double>(efficiency_counters_tracks.track_eta_and_efficiency_counter[layer]) / efficiency_counters_tracks.track_triggered_events_external;

            efficiency_results_tracks.track_eta1_efficiency_external_error[layer] = sqrt(efficiency_results_tracks.track_eta1_efficiency_external[layer] * (1 - efficiency_results_tracks.track_eta1_efficiency_external[layer]) / efficiency_counters_tracks.track_triggered_events_external);
            efficiency_results_tracks.track_eta2_efficiency_external_error[layer] = sqrt(efficiency_results_tracks.track_eta2_efficiency_external[layer] * (1 - efficiency_results_tracks.track_eta2_efficiency_external[layer]) / efficiency_counters_tracks.track_triggered_events_external);
            efficiency_results_tracks.track_eta_or_efficiency_external_error[layer] = sqrt(efficiency_results_tracks.track_eta_or_efficiency_external[layer] * (1 - efficiency_results_tracks.track_eta_or_efficiency_external[layer]) / efficiency_counters_tracks.track_triggered_events_external);
            efficiency_results_tracks.track_eta_and_efficiency_external_error[layer] = sqrt(efficiency_results_tracks.track_eta_and_efficiency_external[layer] * (1 - efficiency_results_tracks.track_eta_and_efficiency_external[layer]) / efficiency_counters_tracks.track_triggered_events_external);
        }
    }
    // External trigger + RPC as a trigger efficiency results
    for (int layer = 0; layer < 3; layer++) {
        if (efficiency_counters.triggered_events_rpc[layer] > 0) {
            efficiency_results.eta1_efficiency_rpc[layer] = static_cast<double>(efficiency_counters.eta1_efficiency_counter_rpc[layer]) / efficiency_counters.triggered_events_rpc[layer];
            efficiency_results.eta2_efficiency_rpc[layer] = static_cast<double>(efficiency_counters.eta2_efficiency_counter_rpc[layer]) / efficiency_counters.triggered_events_rpc[layer];
            efficiency_results.eta_or_efficiency_rpc[layer] = static_cast<double>(efficiency_counters.eta_or_efficiency_counter_rpc[layer]) / efficiency_counters.triggered_events_rpc[layer];
            efficiency_results.eta_and_efficiency_rpc[layer] = static_cast<double>(efficiency_counters.eta_and_efficiency_counter_rpc[layer]) / efficiency_counters.triggered_events_rpc[layer];

            efficiency_results.eta1_efficiency_rpc_error[layer] = sqrt(efficiency_results.eta1_efficiency_rpc[layer] * (1 - efficiency_results.eta1_efficiency_rpc[layer]) / efficiency_counters.triggered_events_rpc[layer]);
            efficiency_results.eta2_efficiency_rpc_error[layer] = sqrt(efficiency_results.eta2_efficiency_rpc[layer] * (1 - efficiency_results.eta2_efficiency_rpc[layer]) / efficiency_counters.triggered_events_rpc[layer]);
            efficiency_results.eta_or_efficiency_rpc_error[layer] = sqrt(efficiency_results.eta_or_efficiency_rpc[layer] * (1 - efficiency_results.eta_or_efficiency_rpc[layer]) / efficiency_counters.triggered_events_rpc[layer]);
            efficiency_results.eta_and_efficiency_rpc_error[layer] = sqrt(efficiency_results.eta_and_efficiency_rpc[layer] * (1 - efficiency_results.eta_and_efficiency_rpc[layer]) / efficiency_counters.triggered_events_rpc[layer]);
        }
        if (efficiency_counters_tracks.track_triggered_events_rpc[layer] > 0) {
            efficiency_results_tracks.track_eta1_efficiency_rpc[layer] = static_cast<double>(efficiency_counters_tracks.track_eta1_efficiency_counter_rpc[layer]) / efficiency_counters_tracks.track_triggered_events_rpc[layer];
            efficiency_results_tracks.track_eta2_efficiency_rpc[layer] = static_cast<double>(efficiency_counters_tracks.track_eta2_efficiency_counter_rpc[layer]) / efficiency_counters_tracks.track_triggered_events_rpc[layer];
            efficiency_results_tracks.track_eta_or_efficiency_rpc[layer] = static_cast<double>(efficiency_counters_tracks.track_eta_or_efficiency_counter_rpc[layer]) / efficiency_counters_tracks.track_triggered_events_rpc[layer];
            efficiency_results_tracks.track_eta_and_efficiency_rpc[layer] = static_cast<double>(efficiency_counters_tracks.track_eta_and_efficiency_counter_rpc[layer]) / efficiency_counters_tracks.track_triggered_events_rpc[layer];

            efficiency_results_tracks.track_eta1_efficiency_rpc_error[layer] = sqrt(efficiency_results_tracks.track_eta1_efficiency_rpc[layer] * (1 - efficiency_results_tracks.track_eta1_efficiency_rpc[layer]) / efficiency_counters_tracks.track_triggered_events_rpc[layer]);
            efficiency_results_tracks.track_eta2_efficiency_rpc_error[layer] = sqrt(efficiency_results_tracks.track_eta2_efficiency_rpc[layer] * (1 - efficiency_results_tracks.track_eta2_efficiency_rpc[layer]) / efficiency_counters_tracks.track_triggered_events_rpc[layer]);
            efficiency_results_tracks.track_eta_or_efficiency_rpc_error[layer] = sqrt(efficiency_results_tracks.track_eta_or_efficiency_rpc[layer] * (1 - efficiency_results_tracks.track_eta_or_efficiency_rpc[layer]) / efficiency_counters_tracks.track_triggered_events_rpc[layer]);
            efficiency_results_tracks.track_eta_and_efficiency_rpc_error[layer] = sqrt(efficiency_results_tracks.track_eta_and_efficiency_rpc[layer] * (1 - efficiency_results_tracks.track_eta_and_efficiency_rpc[layer]) / efficiency_counters_tracks.track_triggered_events_rpc[layer]);
        }
    }
}

/// Utility function to create efficiency histograms in nested directories
void DataAnalyzer::createHistograms() {
    if (!output_file) return;

    // External trigger only efficiency histograms
    output_file->mkdir("efficiencies_histograms/external_trigger");
    output_file->cd("efficiencies_histograms/external_trigger");

    for (int layer = 0; layer < 3; layer++) {
        TH1F* hist_eff_external_trigger = new TH1F(
            Form("eff_external_trigger_layer%d", layer),
            Form("Layer %d Efficiency (External Trigger)", layer),
            4, 0.5, 4.5
        );
        hist_eff_external_trigger->GetXaxis()->SetBinLabel(1, "eta1");
        hist_eff_external_trigger->GetXaxis()->SetBinLabel(2, "eta2");
        hist_eff_external_trigger->GetXaxis()->SetBinLabel(3, "OR");
        hist_eff_external_trigger->GetXaxis()->SetBinLabel(4, "AND");

        hist_eff_external_trigger->SetBinContent(1, efficiency_results.eta1_efficiency_external[layer]);
        hist_eff_external_trigger->SetBinError(1, efficiency_results.eta1_efficiency_external_error[layer]);
        hist_eff_external_trigger->SetBinContent(2, efficiency_results.eta2_efficiency_external[layer]);
        hist_eff_external_trigger->SetBinError(2, efficiency_results.eta2_efficiency_external_error[layer]);
        hist_eff_external_trigger->SetBinContent(3, efficiency_results.eta_or_efficiency_external[layer]);
        hist_eff_external_trigger->SetBinError(3, efficiency_results.eta_or_efficiency_external_error[layer]);
        hist_eff_external_trigger->SetBinContent(4, efficiency_results.eta_and_efficiency_external[layer]);
        hist_eff_external_trigger->SetBinError(4, efficiency_results.eta_and_efficiency_external_error[layer]);

        hist_eff_external_trigger->Write();
    }

    // External trigger + RPC as trigger efficiency histograms
    output_file->cd();
    output_file->mkdir("efficiencies_histograms/external_plus_rpc_trigger");
    output_file->cd("efficiencies_histograms/external_plus_rpc_trigger");
    
    for (int layer = 0; layer < 3; layer++) {
        TH1F* hist_eff_rpc_trigger = new TH1F(
            Form("eff_rpc_layer%d", layer),
            Form("Layer %d Efficiency (External and RPC Trigger)", layer),
            4, 0.5, 4.5
        );
        hist_eff_rpc_trigger->GetXaxis()->SetBinLabel(1, "eta1");
        hist_eff_rpc_trigger->GetXaxis()->SetBinLabel(2, "eta2");
        hist_eff_rpc_trigger->GetXaxis()->SetBinLabel(3, "OR");
        hist_eff_rpc_trigger->GetXaxis()->SetBinLabel(4, "AND");

        hist_eff_rpc_trigger->SetBinContent(1, efficiency_results.eta1_efficiency_rpc[layer]);
        hist_eff_rpc_trigger->SetBinError(1, efficiency_results.eta1_efficiency_rpc_error[layer]);
        hist_eff_rpc_trigger->SetBinContent(2, efficiency_results.eta2_efficiency_rpc[layer]);
        hist_eff_rpc_trigger->SetBinError(2, efficiency_results.eta2_efficiency_rpc_error[layer]);
        hist_eff_rpc_trigger->SetBinContent(3, efficiency_results.eta_or_efficiency_rpc[layer]);
        hist_eff_rpc_trigger->SetBinError(3, efficiency_results.eta_or_efficiency_rpc_error[layer]);
        hist_eff_rpc_trigger->SetBinContent(4, efficiency_results.eta_and_efficiency_rpc[layer]);
        hist_eff_rpc_trigger->SetBinError(4, efficiency_results.eta_and_efficiency_rpc_error[layer]);

        hist_eff_rpc_trigger->Write();
    }

    // Track-based external trigger efficiency histograms
    output_file->cd();
    output_file->mkdir("efficiencies_histograms/track_external_trigger");
    output_file->cd("efficiencies_histograms/track_external_trigger");

    for (int layer = 0; layer < 3; layer++) {
        TH1F* hist_track_eff_external_trigger = new TH1F(
            Form("track_eff_external_trigger_layer%d", layer),
            Form("Layer %d Efficiency (Track and External Trigger)", layer),
            4, 0.5, 4.5
        );
        hist_track_eff_external_trigger->GetXaxis()->SetBinLabel(1, "eta1");
        hist_track_eff_external_trigger->GetXaxis()->SetBinLabel(2, "eta2");
        hist_track_eff_external_trigger->GetXaxis()->SetBinLabel(3, "OR");
        hist_track_eff_external_trigger->GetXaxis()->SetBinLabel(4, "AND");

        hist_track_eff_external_trigger->SetBinContent(1, efficiency_results_tracks.track_eta1_efficiency_external[layer]);
        hist_track_eff_external_trigger->SetBinError(1, efficiency_results_tracks.track_eta1_efficiency_external_error[layer]);
        hist_track_eff_external_trigger->SetBinContent(2, efficiency_results_tracks.track_eta2_efficiency_external[layer]);
        hist_track_eff_external_trigger->SetBinError(2, efficiency_results_tracks.track_eta2_efficiency_external_error[layer]);
        hist_track_eff_external_trigger->SetBinContent(3, efficiency_results_tracks.track_eta_or_efficiency_external[layer]);
        hist_track_eff_external_trigger->SetBinError(3, efficiency_results_tracks.track_eta_or_efficiency_external_error[layer]);
        hist_track_eff_external_trigger->SetBinContent(4, efficiency_results_tracks.track_eta_and_efficiency_external[layer]);
        hist_track_eff_external_trigger->SetBinError(4, efficiency_results_tracks.track_eta_and_efficiency_external_error[layer]);

        hist_track_eff_external_trigger->Write();
    }

    // Track-based external trigger + RPC as trigger efficiency histograms
    output_file->cd();
    output_file->mkdir("efficiencies_histograms/track_external_plus_rpc_trigger");
    output_file->cd("efficiencies_histograms/track_external_plus_rpc_trigger");

    for (int layer = 0; layer < 3; layer++) {
        TH1F* hist_track_eff_rpc_trigger = new TH1F(
            Form("track_eff_rpc_layer%d", layer),
            Form("Layer %d Efficiency (Track and RPC Trigger)", layer),
            4, 0.5, 4.5
        );
        hist_track_eff_rpc_trigger->GetXaxis()->SetBinLabel(1, "eta1");
        hist_track_eff_rpc_trigger->GetXaxis()->SetBinLabel(2, "eta2");
        hist_track_eff_rpc_trigger->GetXaxis()->SetBinLabel(3, "OR");
        hist_track_eff_rpc_trigger->GetXaxis()->SetBinLabel(4, "AND");

        hist_track_eff_rpc_trigger->SetBinContent(1, efficiency_results_tracks.track_eta1_efficiency_rpc[layer]);
        hist_track_eff_rpc_trigger->SetBinError(1, efficiency_results_tracks.track_eta1_efficiency_rpc_error[layer]);
        hist_track_eff_rpc_trigger->SetBinContent(2, efficiency_results_tracks.track_eta2_efficiency_rpc[layer]);
        hist_track_eff_rpc_trigger->SetBinError(2, efficiency_results_tracks.track_eta2_efficiency_rpc_error[layer]);
        hist_track_eff_rpc_trigger->SetBinContent(3, efficiency_results_tracks.track_eta_or_efficiency_rpc[layer]);
        hist_track_eff_rpc_trigger->SetBinError(3, efficiency_results_tracks.track_eta_or_efficiency_rpc_error[layer]);
        hist_track_eff_rpc_trigger->SetBinContent(4, efficiency_results_tracks.track_eta_and_efficiency_rpc[layer]);
        hist_track_eff_rpc_trigger->SetBinError(4, efficiency_results_tracks.track_eta_and_efficiency_rpc_error[layer]);

        hist_track_eff_rpc_trigger->Write();
    }

    // Return to root directory
    output_file->cd();
}


void DataAnalyzer::clearEventVectors() {
    hit_clk.clear();
    hit_channel.clear();
    hit_raw_bcid.clear();
    hit_bcid.clear();
    hit_time1.clear();
    hit_time2.clear();
    hit_rise.clear();
    hit_raw_bcout.clear();
    hit_bcout.clear();

    proc_layer.clear();
    proc_strip.clear();
    proc_time1.clear();
    proc_time2.clear();
    proc_dt_time1_time2.clear();
    proc_trigger_time.clear();
    proc_dt_time1_trigger.clear();
    proc_dt_time2_trigger.clear();
    proc_tot1.clear();
    proc_tot2.clear();

    cluster_size_eta1.clear();
    cluster_size_eta2.clear();
    cluster_tot1.clear();
    cluster_tot2.clear();

    track_length_eta1.clear();
    track_length_eta2.clear();
    track_width_eta1.clear();
    track_width_eta2.clear();
    track_size_eta1.clear();
    track_size_eta2.clear();
}

/// Utility function to push raw hit data into the corresponding vectors for tree filling
void DataAnalyzer::pushBackHitData(const Hit& hit) {
    hit_clk.push_back(hit.getClk());
    hit_channel.push_back(hit.getChannel());
    hit_raw_bcid.push_back(hit.getRawBCID());
    hit_bcid.push_back(hit.getBCID());
    hit_time1.push_back(hit.getRawTimeEta1());
    hit_time2.push_back(hit.getRawTimeEta2());
    hit_rise.push_back(hit.getRise());
    hit_raw_bcout.push_back(hit.getRawBCOut());
    hit_bcout.push_back(hit.getBCOut());
}

/// Utility function to push processed hit data into the corresponding vectors for tree filling
void DataAnalyzer::pushBackProcessedData(const Event& event) {

    proc_trigger_time.push_back(event.getTriggerTime());
    // Loop over hits in the event and push processed data into vectors
    for (const auto& hit : event.getHits()) {
        proc_layer.push_back(hit.getLayer());
        proc_strip.push_back(hit.getStrip());
        proc_time1.push_back(hit.getTimeEta1());
        proc_time2.push_back(hit.getTimeEta2());
        proc_dt_time1_time2.push_back(hit.getTimeEta1() - hit.getTimeEta2());
        proc_dt_time1_trigger.push_back(hit.getTimeEta1() - event.getTriggerTime());
        proc_dt_time2_trigger.push_back(hit.getTimeEta2() - event.getTriggerTime());
        proc_tot1.push_back(hit.getToT1());
        proc_tot2.push_back(hit.getToT2());
    }
    /// WIP: Decide whether to do filtering first: 
    /// FILTERING: Only push information of the non-trigger channel rising hits with valid time information (to match with ToT calculation which only processes rising edges)
    /*
    for (const auto& hit : event.getHits()) {
        // Only process non-trigger channel and rising edges for consistency with ToT calculation
        if (hit.getChannel() == event.getTriggerChannel() || hit.getRise() != 1) continue;
        
        if (hit.hasEta1Time() || hit.hasEta2Time()) {
            proc_layer.push_back(hit.getLayer());
            proc_strip.push_back(hit.getStrip());
        }
        if (hit.hasEta1Time()) {
            proc_time1.push_back(hit.getTimeEta1());
            proc_dt_time1_trigger.push_back(hit.getTimeEta1() - event.getTriggerTime());
        }
        if (hit.hasEta2Time()) {
            proc_time2.push_back(hit.getTimeEta2());
            proc_dt_time2_trigger.push_back(hit.getTimeEta2() - event.getTriggerTime());
        }
        if (hit.hasEta1Time() && hit.hasEta2Time()) proc_dt_time1_time2.push_back(hit.getTimeEta1() - hit.getTimeEta2());
        if (hit.getToT1() > 0) proc_tot1.push_back(hit.getToT1());
        if (hit.getToT2() > 0) proc_tot2.push_back(hit.getToT2());
    }
    */
}

/// Utility functions to push cluster-level data into the corresponding vectors for tree filling
void DataAnalyzer::pushBackClusterDataEta1(const Cluster& cluster) {
    cluster_size_eta1.push_back(cluster.getSize()); // Total cluster size
    if (cluster.getTot1() > 0) cluster_tot1.push_back(cluster.getTot1());
}

void DataAnalyzer::pushBackClusterDataEta2(const Cluster& cluster) {
    cluster_size_eta2.push_back(cluster.getSize()); // Total cluster size
    if (cluster.getTot2() > 0) cluster_tot2.push_back(cluster.getTot2());
}

/// Utility functions to push track-level data into the corresponding vectors for tree filling
void DataAnalyzer::pushBackTrackDataEta1(const Track& track) {
    track_length_eta1.push_back(track.getLayerCount(0)); // Eta1 layer count
    track_width_eta1.push_back(track.getWidth(0));       // Eta1 width
    track_size_eta1.push_back(track.getSize(0));         // Eta1 size
}

void DataAnalyzer::pushBackTrackDataEta2(const Track& track) {
    track_length_eta2.push_back(track.getLayerCount(1)); // Eta2 layer count
    track_width_eta2.push_back(track.getWidth(1));       // Eta2 width
    track_size_eta2.push_back(track.getSize(1));         // Eta2 size
}

/// Main entry point for processing input data from the input file
void DataAnalyzer::processInputData(const std::string& file_path, const int dt_max_arg, const int dt_min_arg) {
    // Store time window parameters for use throughout processing
    dt_max = dt_max_arg;
    dt_min = dt_min_arg;
    
    namespace fs = std::filesystem;

    if (!fs::exists(file_path)) {
        std::cerr << "ERROR: Input path does not exist: " << file_path << std::endl;
        return;
    }

    std::cout << "Processing file: " << file_path << std::endl;
    processFile(file_path);
}

/// Process a single file by reading lines and extracting word data
void DataAnalyzer::processFile(const std::string& file_path) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open file: " << file_path << std::endl;
        return;
    }

    int clk, word, raw_bcout;

    // Initialize efficiency counters and struct once before processing all words in this file
    // Use member variables so they persist after processing
    efficiency_counters = {};
    efficiency_counters_tracks = {};

    // Read words from file: "clk word raw_bcout" format (all in hexadecimal or as they appear)
    while (infile >> std::dec >> clk >> std::hex >> word >> raw_bcout) {
        processSingleWord(clk, word, raw_bcout, efficiency_counters, efficiency_counters_tracks);
    }

    // Process any remaining event at the end of file (rare, but handles incomplete final events)
    if (current_event_hits.size() > 0) {
        processEvent(efficiency_counters, efficiency_counters_tracks);
        current_event_number++;
        current_event_hits.clear();
        clearEventVectors();
    }

    // Update efficiencies from counters before writing histograms
    updateEfficiencies();
    
    // Create and write efficiency histograms to output file
    createHistograms();

    infile.close();
}

/// Helper function to extract raw BCID from DCT word
int DataAnalyzer::extractRawBCID(int word) {
    int rise = word & 0x01;               // Bit 0
    int raw_bcid;
    
    if (rise == 1) {  // Rising edge
        raw_bcid = (word >> 12) & 0xFF;   // Bits 12-19
    } else {          // Falling edge
        raw_bcid = (word >> 11) & 0x1FF;  // Bits 11-19
    }
    
    return raw_bcid;
}

/// Process a single word and accumulate into events
/// When CLK reaches 127 (end of BC frame), process the accumulated event
void DataAnalyzer::processSingleWord(int clk, int word, int raw_bcout, EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks) {
    // Only process not empty words
    if (word != EMPTY_WORD) {
        // Calculate BC0 if this is the first hit of an event
        if (current_event_hits.empty()) {
            int raw_bcid = extractRawBCID(word);
            BC0 = raw_bcid % 256;
        }

        // Create a Hit object from the raw word data
        Hit hit(clk, word, raw_bcout, BC0);
        hit.setIdx(current_event_hits.size());
        current_event_hits.push_back(hit);

        // Accumulate raw hit data into vector
        pushBackHitData(hit);
    }

    // Check for end of event and process it
    if (clk == 127 && current_event_hits.size() > 0) {
        processEvent(counters, counters_tracks);
        current_event_number++;
        current_event_hits.clear();

        // Clear processed data vectors for the next event
        clearEventVectors();
    }
}

/// Process a complete event that has been accumulated
void DataAnalyzer::processEvent(EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks) {
    n_hits = current_event_hits.size();
    std::cout << "\n" << std::endl;
    std::cout << "╔═════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║ Processing event No. " << current_event_number << " with " << n_hits << " hits" << std::endl;
    
    // Fill the input data tree with raw hit information
    if (input_data_tree && hit_clk.size() > 0) {
        input_data_tree->Fill();
    }

    // Create an Event object and transfer hits using move semantics, passing the counters reference
    Event event(current_event_number, std::move(current_event_hits), counters, counters_tracks);

    // Extract trigger information from the event
    event.extractTriggerTime();

    // Calculate ToT before clusterization
    event.calculateTOT();

    // Push back processed hit data into vectors for tree filling
    pushBackProcessedData(event);

    // Fill the processed data tree
    if (processed_data_tree && proc_layer.size() > 0) {
        processed_data_tree->Fill();
    }

    // Clusterization
    event.clusterize();

    // Calculate ToT for cluster centers
    event.calculateTOTCluster();

    // Push back cluster-level data for both eta1 and eta2 sides
    for (const auto& cluster : event.getClustersEta1()) {
        pushBackClusterDataEta1(cluster);
    }
    for (const auto& cluster : event.getClustersEta2()) {
        pushBackClusterDataEta2(cluster);
    }

    // Fill the clusterization tree (write if either eta1 or eta2 has clusters)
    if (clusterization_tree && (cluster_size_eta1.size() > 0 && cluster_size_eta2.size() > 0)) {
        clusterization_tree->Fill();
    }

    // Track reconstruction
    event.reconstructTracks();

    // Push back track-level data for both eta1 and eta2 sides
    for (const auto& track : event.getTracksEta1()) {
        pushBackTrackDataEta1(track);
    }
    for (const auto& track : event.getTracksEta2()) {
        pushBackTrackDataEta2(track);
    }

    // Reset and update efficiency flags for each new event
    event.updateEfficiencyFlags(dt_max, dt_min);

    // Update efficiency counters based on reconstructed tracks
    event.updateEfficiencyCounters();

    // Fill the track reconstruction tree (write if either eta1 or eta2 has tracks)
    if (track_reconstruction_tree && (track_length_eta1.size() > 0 && track_length_eta2.size() > 0)) {
        track_reconstruction_tree->Fill();
    }
}