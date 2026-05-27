#include "dataProcesser.hpp"


/// TODO:
/// - Channel inversion for eta2 (?) side (in the same column: channel 1 of swaps with 8, 2 with 7, etc.)
/// - Lift restrictions on track length/size to have them be of different sizes

// ==========================================================================================
// DataProcesser class implementation for processing raw DCT data
// ==========================================================================================
/// Constructor and destructor for DataProcesser class
DataProcesser::DataProcesser() {
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
    use_external_trigger = true;
    initializeCounters();
}

DataProcesser::~DataProcesser() {
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
void DataProcesser::setupOutputFile() {
    output_file = new TFile("output.root", "RECREATE");
    input_data_tree = new TTree("InputData", "Raw hit data from the DCT");
    processed_data_tree = new TTree("ProcessedData", "Processed hit and event-level data");
    clusterization_tree = new TTree("Clusterization", "Cluster-level data");
    track_reconstruction_tree = new TTree("TrackReconstruction", "Track-level data");
}

/// Utility function to setup branches for all trees in the output file
void DataProcesser::setupBranches() {
    // Branch definitions for raw data tree
    input_data_tree->Branch("hit_clk", &hit_clk);
    input_data_tree->Branch("hit_channel", &hit_channel);
    input_data_tree->Branch("hit_raw_bcid", &hit_raw_bcid);
    input_data_tree->Branch("bc0", &BC0);
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
    clusterization_tree->Branch("cluster_size_eta1_layer0", &cluster_size_eta1_layers[0]);
    clusterization_tree->Branch("cluster_size_eta1_layer1", &cluster_size_eta1_layers[1]);
    clusterization_tree->Branch("cluster_size_eta1_layer2", &cluster_size_eta1_layers[2]);
    clusterization_tree->Branch("cluster_size_eta2_layer0", &cluster_size_eta2_layers[0]);
    clusterization_tree->Branch("cluster_size_eta2_layer1", &cluster_size_eta2_layers[1]);
    clusterization_tree->Branch("cluster_size_eta2_layer2", &cluster_size_eta2_layers[2]);
    clusterization_tree->Branch("cluster_tot1_layer0", &cluster_tot1_layers[0]);
    clusterization_tree->Branch("cluster_tot1_layer1", &cluster_tot1_layers[1]);
    clusterization_tree->Branch("cluster_tot1_layer2", &cluster_tot1_layers[2]);
    clusterization_tree->Branch("cluster_tot2_layer0", &cluster_tot2_layers[0]);
    clusterization_tree->Branch("cluster_tot2_layer1", &cluster_tot2_layers[1]);
    clusterization_tree->Branch("cluster_tot2_layer2", &cluster_tot2_layers[2]);

    // Branch definitions for track reconstruction tree
    track_reconstruction_tree->Branch("track_length_eta1", &track_length_eta1);
    track_reconstruction_tree->Branch("track_length_eta2", &track_length_eta2);
    track_reconstruction_tree->Branch("track_width_eta1", &track_width_eta1);
    track_reconstruction_tree->Branch("track_width_eta2", &track_width_eta2);
    track_reconstruction_tree->Branch("track_size_eta1", &track_size_eta1);
    track_reconstruction_tree->Branch("track_size_eta2", &track_size_eta2);
}

/// Utility function to initialize efficiency counters and results
void DataProcesser::initializeCounters() {
    efficiency_counters = {};
    efficiency_counters_tracks = {};
    efficiency_results = {};
    efficiency_results_tracks = {};
}

/// Efficiency live calculation based on reconstructed tracks
void DataProcesser::updateEfficiencies() {
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
void DataProcesser::createHistograms() {
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

/// Utility function to clear all vectors for the next event processing
void DataProcesser::clearEventVectors() {
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
    for (int layer = 0; layer < 3; layer++) {
        cluster_size_eta1_layers[layer].clear();
        cluster_size_eta2_layers[layer].clear();
        cluster_tot1_layers[layer].clear();
        cluster_tot2_layers[layer].clear();
    }

    track_length_eta1.clear();
    track_length_eta2.clear();
    track_width_eta1.clear();
    track_width_eta2.clear();
    track_size_eta1.clear();
    track_size_eta2.clear();
}

/// Utility function to push raw hit data into the corresponding vectors for tree filling
void DataProcesser::pushBackHitData(const Hit& hit) {
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
void DataProcesser::pushBackProcessedData(const Event& event) {
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

    /*
    /// WIP: Decide whether to do filtering first: 
    /// FILTERING: Only push information of the non-trigger channel rising hits with valid time information (to match with ToT calculation which only processes rising edges)
    for (const auto& hit : event.getHits()) {
        // Only process non-trigger channel and rising edges for consistency with ToT calculation
        if (hit.getChannel() == event.getTriggerChannel() || hit.getRise() != 1) continue;
        
        if (hit.hasEta1Time() || hit.hasEta2Time()) {
            proc_layer.push_back(hit.getLayer());
            proc_strip.push_back(hit.getStrip());
        }
        if (hit.hasEta1Time()) {
            if (hit.getTimeEta1() < -1) {
                // std::cout << "FATAL: Hit number " << hit.getIdx() << " with negative BCID after BC0 correction: " << hit.getBCID() << " (raw BCID: " << hit.getRawBCID() << ", BC0: " << BC0 << ")" << std::endl;
                // std::exit(1);
            }
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
void DataProcesser::pushBackClusterData(const Cluster& cluster) {
    if (cluster.ETA1) {
        // Cluster size information
        cluster_size_eta1.push_back(cluster.getSize());
        cluster_size_eta1_layers[cluster.getLayer()].push_back(cluster.getSize());

        // ToT information
        if (cluster.getTot1() > 0) {
            cluster_tot1.push_back(cluster.getTot1());
            cluster_tot1_layers[cluster.getLayer()].push_back(cluster.getTot1());
        }
    } else if (cluster.ETA2) { // ETA2 side, same structure but different vectors
        cluster_size_eta2.push_back(cluster.getSize());
        cluster_size_eta2_layers[cluster.getLayer()].push_back(cluster.getSize());
        if (cluster.getTot2() > 0) {
            cluster_tot2.push_back(cluster.getTot2());
            cluster_tot2_layers[cluster.getLayer()].push_back(cluster.getTot2());
        }
    }
}


/// Utility functions to push track-level data into the corresponding vectors for tree filling
void DataProcesser::pushBackTrackDataEta1(const Track& track) {
    track_length_eta1.push_back(track.getLayerCount(0)); // Eta1 layer count
    track_width_eta1.push_back(track.getWidth(0));       // Eta1 width
    track_size_eta1.push_back(track.getSize(0));         // Eta1 size
}

void DataProcesser::pushBackTrackDataEta2(const Track& track) {
    track_length_eta2.push_back(track.getLayerCount(1)); // Eta2 layer count
    track_width_eta2.push_back(track.getWidth(1));       // Eta2 width
    track_size_eta2.push_back(track.getSize(1));         // Eta2 size
}

/// Main entry point for processing input data from the input file
void DataProcesser::processInputData(const std::string& file_path, const int dt_max_arg, const int dt_min_arg, InputFormat format, bool use_external_trigger_arg) {
    // Store time window parameters for use throughout processing
    dt_max = dt_max_arg;
    dt_min = dt_min_arg;
    use_external_trigger = use_external_trigger_arg;
    
    namespace fs = std::filesystem;

    if (!fs::exists(file_path)) {
        std::cerr << "ERROR: Input path does not exist: " << file_path << std::endl;
        return;
    }

    std::cout << "Processing file: " << file_path << std::endl;
    if (format == InputFormat::DecodedWords) {
        processFileDecoded(file_path);
    } else {
        processFileFiledump(file_path);
    }
}

/// Process a single file by reading lines and extracting word data
/// New format: filedump_*.txt packets as produced by fast DAQ
void DataProcesser::processFileFiledump(const std::string& file_path) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open file: " << file_path << std::endl;
        return;
    }

    const int raw_bcout = 0;
    int clk = 0;
    bool processing_event = false;
    std::string token;
    std::vector<int> event_words;
    std::vector<int> event_clks;
    int min_bcid_event = 256;

    // Initialize efficiency counters and struct once before processing all words in this file
    // Use member variables so they persist after processing
    efficiency_counters = {};
    efficiency_counters_tracks = {};

    while (infile >> token) {
        if (token.length() <= 512) {
            continue;
        }

        std::string header = token.substr(0, 2);
        if (header == "0a") {
            // PASS 2: Process all words from this event using event-level BC0
            if (processing_event && !event_words.empty()) {
                if (min_bcid_event < 256) {
                    BC0 = min_bcid_event;
                }
                for (size_t i = 0; i < event_words.size(); ++i) {
                    processSingleWord(event_clks[i], event_words[i], raw_bcout, efficiency_counters, efficiency_counters_tracks, false);
                }
                if (current_event_hits.size() > 0) {
                    processEvent(efficiency_counters, efficiency_counters_tracks);
                    current_event_number++;
                    current_event_hits.clear();
                    clearEventVectors();
                }
            }

            processing_event = true;
            clk = 0;
            event_words.clear();
            event_clks.clear();
            min_bcid_event = 256;
        }

        if (!processing_event) {
            continue;
        }

        // PASS 1: Scan all words in the event to find the minimum BCID (BC0)
        for (int i = 0; i < 2 * 128; i++) {
            const size_t offset = static_cast<size_t>(i) * 8 + 2;
            if (offset + 8 > token.length()) {
                break;
            }

            unsigned int word = 0;
            bool valid_hex = true;
            for (int j = 0; j < 8; ++j) {
                char c = token[offset + j];
                word <<= 4;
                if (c >= '0' && c <= '9') {
                    word |= static_cast<unsigned int>(c - '0');
                } else if (c >= 'a' && c <= 'f') {
                    word |= static_cast<unsigned int>(c - 'a' + 10);
                } else if (c >= 'A' && c <= 'F') {
                    word |= static_cast<unsigned int>(c - 'A' + 10);
                } else {
                    valid_hex = false;
                    break;
                }
            }

            if (!valid_hex) {
                continue;
            }

            int dct = (word >> 28) & 0xF;
            if (dct == 1) {
                clk++;
            }

            if ((word & 0x0FFFFFFF) == 0x05555555) {
                continue;
            }

            event_clks.push_back(clk);
            event_words.push_back(static_cast<int>(word));

            int raw_bcid = extractRawBCID(static_cast<int>(word));
            int raw_bcid_mod = raw_bcid & 0xFF;
            if (raw_bcid_mod < min_bcid_event) {
                min_bcid_event = raw_bcid_mod;
            }
        }
    }

    // Process any remaining event at the end of file
    if (!event_words.empty()) {
        if (min_bcid_event < 256) {
            BC0 = min_bcid_event;
        }
        // PASS 2: Process all words from this event using event-level BC0
        for (size_t i = 0; i < event_words.size(); ++i) {
            processSingleWord(event_clks[i], event_words[i], raw_bcout, efficiency_counters, efficiency_counters_tracks, false);
        }
        if (current_event_hits.size() > 0) {
            processEvent(efficiency_counters, efficiency_counters_tracks);
            current_event_number++;
            current_event_hits.clear();
            clearEventVectors();
        }
    }

    // Update efficiencies from counters before writing histograms
    updateEfficiencies();

    // Create and write efficiency histograms to output file
    createHistograms();

    infile.close();
}

/// Process a single file by reading decoded words: "clk word raw_bcout" format
void DataProcesser::processFileDecoded(const std::string& file_path) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open file: " << file_path << std::endl;
        return;
    }

    int clk = 0;
    int word = 0;
    int raw_bcout = 0;

    // Initialize efficiency counters and struct once before processing all words in this file
    // Use member variables so they persist after processing
    efficiency_counters = {};
    efficiency_counters_tracks = {};

    // Read words from file: "clk word raw_bcout" format
    while (infile >> std::dec >> clk >> std::hex >> word >> raw_bcout) {
        processSingleWord(clk, word, raw_bcout, efficiency_counters, efficiency_counters_tracks);
    }

    // Process any remaining event at the end of file
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
int DataProcesser::extractRawBCID(int word) {
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
void DataProcesser::processSingleWord(int clk, int word, int raw_bcout, EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks, bool end_on_clk) {
    // Only process not empty words
    if (word != EMPTY_WORD) {
        // Create a Hit object from the raw word data
        Hit hit(clk, word, raw_bcout, BC0);
        hit.setIdx(current_event_hits.size());
        current_event_hits.push_back(hit);

        // Accumulate raw hit data into vector
        pushBackHitData(hit);
    }

    // Check for end of event and process it
    if (end_on_clk && clk == 127 && current_event_hits.size() > 0) {
        processEvent(counters, counters_tracks);
        current_event_number++;
        current_event_hits.clear();

        // Clear processed data vectors for the next event
        clearEventVectors();
    }
}

/// Process a complete event that has been accumulated
void DataProcesser::processEvent(EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks) {
    n_hits = current_event_hits.size();
    // std::cout << "\n" << std::endl;
    // std::cout << "╔═════════════════════════════════════════════════════════════════╗" << std::endl;
    // std::cout << "║ Processing event No. " << current_event_number << " with " << n_hits << " hits" << std::endl;
    
    // Fill the input data tree with raw hit information
    if (input_data_tree && hit_clk.size() > 0) {
        input_data_tree->Fill();
    }

    // Create an Event object and transfer hits using move semantics, passing the counters reference
    Event event(current_event_number, std::move(current_event_hits), counters, counters_tracks, use_external_trigger);

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
        pushBackClusterData(cluster);
    }
    for (const auto& cluster : event.getClustersEta2()) {
        pushBackClusterData(cluster);
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