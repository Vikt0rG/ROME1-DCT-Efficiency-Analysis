#include "dataProcesser.hpp"


/// TODO:
/// - Channel inversion for eta2 (?) side (in the same column: channel 1 of swaps with 8, 2 with 7, etc.)
/// - Lift restrictions on track length/size to have them be of different sizes

// ==========================================================================================
// DataProcesser class implementation for processing raw DCT data
// ==========================================================================================
// Constructor and destructor for DataProcesser class
DataProcesser::DataProcesser(const std::string& input_path, const int dt_max, const int dt_min, InputFormat format, bool use_external_trigger_arg, bool reject_background_arg) 
    : _input_path(input_path)
    , _dt_max(dt_max)
    , _dt_min(dt_min)
    , _format(format)
    , _use_external_trigger(use_external_trigger_arg)
    , _reject_background(reject_background_arg) {
    setupOutputFile();
    setupBranches();
}

DataProcesser::~DataProcesser() {
    if (_output_file) { 
        // Write all trees to file before closing
        if (_input_data_tree) _input_data_tree->Write();
        if (_background_rejection_tree) _background_rejection_tree->Write();
        if (_processed_data_tree) _processed_data_tree->Write();
        if (_clusterization_tree) _clusterization_tree->Write();
        if (_track_reconstruction_tree) _track_reconstruction_tree->Write();
        _output_file->Close(); 
        delete _output_file; 
    }
}

// ------------------------------------------------------------------------------------------
// File, tree and branch setup utility functions

// Utility function to setup the output ROOT file and create the necessary trees
void DataProcesser::setupOutputFile() {
    _output_file = new TFile("output.root", "RECREATE");
    _input_data_tree = new TTree("InputData", "Raw hit data from the DCT", 1, _output_file);
    _processed_data_tree = new TTree("ProcessedData", "Processed hit and event-level data", 1, _output_file);
    _clusterization_tree = new TTree("Clusterization", "Cluster-level data", 1, _output_file);
    _track_reconstruction_tree = new TTree("TrackReconstruction", "Track-level data", 1, _output_file);
}

// Utility function to setup branches for all trees in the output file
void DataProcesser::setupBranches() {
    // Branch definitions for raw data tree
    _input_data_tree->Branch("hit_clk", &hit_clk);
    _input_data_tree->Branch("hit_channel", &hit_channel);
    _input_data_tree->Branch("hit_raw_bcid", &hit_raw_bcid);
    _input_data_tree->Branch("hit_raw_time1", &hit_raw_time1);
    _input_data_tree->Branch("hit_raw_time2", &hit_raw_time2);
    _input_data_tree->Branch("hit_rise", &hit_rise);

    // Branch definitions for processed data tree
    _processed_data_tree->Branch("n_events", &current_event_number);
    _processed_data_tree->Branch("n_hits", &n_hits);
    _processed_data_tree->Branch("proc_layer", &proc_layer);
    _processed_data_tree->Branch("proc_strip", &proc_strip);
    _processed_data_tree->Branch("proc_bc0", &proc_bc0);
    _processed_data_tree->Branch("proc_bcid", &proc_bcid);
    _processed_data_tree->Branch("proc_time1", &proc_time1);
    _processed_data_tree->Branch("proc_time2", &proc_time2);
    _processed_data_tree->Branch("proc_dt_time1_time2", &proc_dt_time1_time2);
    _processed_data_tree->Branch("proc_trigger_time", &proc_trigger_time);
    _processed_data_tree->Branch("proc_dt_time1_trigger", &proc_dt_time1_trigger);
    _processed_data_tree->Branch("proc_dt_time2_trigger", &proc_dt_time2_trigger);
    _processed_data_tree->Branch("proc_tot1", &proc_tot1);
    _processed_data_tree->Branch("proc_tot2", &proc_tot2);

    // Branch definitions for clusterization tree
    _clusterization_tree->Branch("cluster_id_from_eta1", &_cluster_id_from_eta1);
    _clusterization_tree->Branch("cluster_id_from_eta2", &_cluster_id_from_eta2);
    _clusterization_tree->Branch("cluster_size_eta1", &cluster_size_eta1);
    _clusterization_tree->Branch("cluster_size_eta2", &cluster_size_eta2);
    _clusterization_tree->Branch("cluster_tot1_from_eta1", &cluster_tot1_from_eta1);
    _clusterization_tree->Branch("cluster_tot2_from_eta1", &cluster_tot2_from_eta1);
    _clusterization_tree->Branch("cluster_tot1_from_eta2", &cluster_tot1_from_eta2);
    _clusterization_tree->Branch("cluster_tot2_from_eta2", &cluster_tot2_from_eta2);
    _clusterization_tree->Branch("cluster_size_eta1_layer0", &cluster_size_eta1_layers[0]);
    _clusterization_tree->Branch("cluster_size_eta1_layer1", &cluster_size_eta1_layers[1]);
    _clusterization_tree->Branch("cluster_size_eta1_layer2", &cluster_size_eta1_layers[2]);
    _clusterization_tree->Branch("cluster_size_eta2_layer0", &cluster_size_eta2_layers[0]);
    _clusterization_tree->Branch("cluster_size_eta2_layer1", &cluster_size_eta2_layers[1]);
    _clusterization_tree->Branch("cluster_size_eta2_layer2", &cluster_size_eta2_layers[2]);
    _clusterization_tree->Branch("cluster_tot1_from_eta1_layer0", &cluster_tot1_from_eta1_layers[0]);
    _clusterization_tree->Branch("cluster_tot1_from_eta1_layer1", &cluster_tot1_from_eta1_layers[1]);
    _clusterization_tree->Branch("cluster_tot1_from_eta1_layer2", &cluster_tot1_from_eta1_layers[2]);
    _clusterization_tree->Branch("cluster_tot2_from_eta1_layer0", &cluster_tot2_from_eta1_layers[0]);
    _clusterization_tree->Branch("cluster_tot2_from_eta1_layer1", &cluster_tot2_from_eta1_layers[1]);
    _clusterization_tree->Branch("cluster_tot2_from_eta1_layer2", &cluster_tot2_from_eta1_layers[2]);
    _clusterization_tree->Branch("cluster_tot1_from_eta2_layer0", &cluster_tot1_from_eta2_layers[0]);
    _clusterization_tree->Branch("cluster_tot1_from_eta2_layer1", &cluster_tot1_from_eta2_layers[1]);
    _clusterization_tree->Branch("cluster_tot1_from_eta2_layer2", &cluster_tot1_from_eta2_layers[2]);
    _clusterization_tree->Branch("cluster_tot2_from_eta2_layer0", &cluster_tot2_from_eta2_layers[0]);
    _clusterization_tree->Branch("cluster_tot2_from_eta2_layer1", &cluster_tot2_from_eta2_layers[1]);
    _clusterization_tree->Branch("cluster_tot2_from_eta2_layer2", &cluster_tot2_from_eta2_layers[2]);

    // Branch definitions for track reconstruction tree
    _track_reconstruction_tree->Branch("track_id_from_eta1", &_track_id_from_eta1);
    _track_reconstruction_tree->Branch("track_id_from_eta2", &_track_id_from_eta2);
    _track_reconstruction_tree->Branch("in_valid_track_eta1", &_in_valid_track_eta1);
    _track_reconstruction_tree->Branch("in_valid_track_eta2", &_in_valid_track_eta2);
    _track_reconstruction_tree->Branch("track_length_eta1", &track_length_eta1);
    _track_reconstruction_tree->Branch("track_length_eta2", &track_length_eta2);
    _track_reconstruction_tree->Branch("track_width_eta1", &track_width_eta1);
    _track_reconstruction_tree->Branch("track_width_eta2", &track_width_eta2);
    _track_reconstruction_tree->Branch("track_size_eta1", &track_size_eta1);
    _track_reconstruction_tree->Branch("track_size_eta2", &track_size_eta2);
}

// ------------------------------------------------------------------------------------------
// Efficiency calculation helper functions

// Efficiency live calculation based on reconstructed tracks
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

// ------------------------------------------------------------------------------------------
// Histogram creation utility function

// Utility function to create efficiency histograms in nested directories
void DataProcesser::createHistograms() {
    if (!_output_file) return;

    // External trigger only efficiency histograms
    _output_file->mkdir("efficiencies_histograms/external_trigger");
    _output_file->cd("efficiencies_histograms/external_trigger");

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
    _output_file->cd();
    _output_file->mkdir("efficiencies_histograms/external_plus_rpc_trigger");
    _output_file->cd("efficiencies_histograms/external_plus_rpc_trigger");
    
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
    _output_file->cd();
    _output_file->mkdir("efficiencies_histograms/track_external_trigger");
    _output_file->cd("efficiencies_histograms/track_external_trigger");

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
    _output_file->cd();
    _output_file->mkdir("efficiencies_histograms/track_external_plus_rpc_trigger");
    _output_file->cd("efficiencies_histograms/track_external_plus_rpc_trigger");

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
    _output_file->cd();
}

// ------------------------------------------------------------------------------------------
// Clearing and pushing back data utility functions

// Utility function to clear raw data vectors
void DataProcesser::clearRawDataVectors() {
    hit_clk.clear();
    hit_channel.clear();
    hit_raw_bcid.clear();
    hit_raw_time1.clear();
    hit_raw_time2.clear();
    hit_rise.clear();
}

// Utility function to clear all vectors for the next event processing
void DataProcesser::clearEventVectors() {
    proc_layer.clear();
    proc_strip.clear();
    proc_bc0.clear();
    proc_bcid.clear();
    proc_time1.clear();
    proc_time2.clear();
    proc_dt_time1_time2.clear();
    proc_trigger_time.clear();
    proc_dt_time1_trigger.clear();
    proc_dt_time2_trigger.clear();
    proc_tot1.clear();
    proc_tot2.clear();

    _cluster_id_from_eta1.clear();
    _cluster_id_from_eta2.clear();
    cluster_size_eta1.clear();
    cluster_size_eta2.clear();
    cluster_tot1_from_eta1.clear();
    cluster_tot2_from_eta1.clear();
    cluster_tot1_from_eta2.clear();
    cluster_tot2_from_eta2.clear();
    for (int layer = 0; layer < 3; layer++) {
        cluster_size_eta1_layers[layer].clear();
        cluster_size_eta2_layers[layer].clear();
        cluster_tot1_from_eta1_layers[layer].clear();
        cluster_tot2_from_eta1_layers[layer].clear();
        cluster_tot1_from_eta2_layers[layer].clear();
        cluster_tot2_from_eta2_layers[layer].clear();
    }

    _track_id_from_eta1.clear();
    _track_id_from_eta2.clear();
    _in_valid_track_eta1.clear();
    _in_valid_track_eta2.clear();
    track_length_eta1.clear();
    track_length_eta2.clear();
    track_width_eta1.clear();
    track_width_eta2.clear();
    track_size_eta1.clear();
    track_size_eta2.clear();
}

// Utility function to push back raw hit data
void DataProcesser::pushBackWordData(const DCTWord& word) {
    hit_clk.push_back(word.clk);
    hit_channel.push_back(word.channel);
    hit_raw_bcid.push_back(word.raw_bcid);
    hit_raw_time1.push_back(word.raw_time_eta1);
    hit_raw_time2.push_back(word.raw_time_eta2);
    hit_rise.push_back(word.rise);
}

// Utility function to push processed hit data into the corresponding vectors for tree filling
void DataProcesser::pushBackProcessedData(const Event& event) {
    proc_trigger_time.push_back(event.getTriggerTime());
    proc_bc0.push_back(BC0);
    // Loop over hits in the event and push processed data into vectors
    for (const auto& hit : event.getHits()) {
        proc_layer.push_back(hit.getLayer());
        proc_strip.push_back(hit.getStrip());
        proc_bcid.push_back(hit.getBCID());
        proc_time1.push_back(hit.getTimeEta1());
        proc_time2.push_back(hit.getTimeEta2());
        proc_dt_time1_time2.push_back(hit.getTimeEta1() - hit.getTimeEta2());
        proc_dt_time1_trigger.push_back(hit.getTimeEta1() - event.getTriggerTime());
        proc_dt_time2_trigger.push_back(hit.getTimeEta2() - event.getTriggerTime());
        proc_tot1.push_back(hit.getToT1());
        proc_tot2.push_back(hit.getToT2());
    }
}

// Utility functions to push cluster-level data into the corresponding vectors for tree filling
void DataProcesser::pushBackClusterData(const Cluster& cluster) {
    if (cluster.getSide() == Cluster::ETA1) {
        // Cluster size information
        cluster_size_eta1.push_back(cluster.getSize());
        cluster_size_eta1_layers[cluster.getLayer()].push_back(cluster.getSize());

        // ToT information
        if (cluster.getTot1() > 0) {
            cluster_tot1_from_eta1.push_back(cluster.getTot1());
            cluster_tot1_from_eta1_layers[cluster.getLayer()].push_back(cluster.getTot1());
        }
        if (cluster.getTot2() > 0) {
            cluster_tot2_from_eta1.push_back(cluster.getTot2());
            cluster_tot2_from_eta1_layers[cluster.getLayer()].push_back(cluster.getTot2());
        }
    } else if (cluster.getSide() == Cluster::ETA2) { // ETA2 side, same structure but different vectors
        cluster_size_eta2.push_back(cluster.getSize());
        cluster_size_eta2_layers[cluster.getLayer()].push_back(cluster.getSize());
        if (cluster.getTot1() > 0) {
            cluster_tot1_from_eta2.push_back(cluster.getTot1());
            cluster_tot1_from_eta2_layers[cluster.getLayer()].push_back(cluster.getTot1());
        }
        if (cluster.getTot2() > 0) {
            cluster_tot2_from_eta2.push_back(cluster.getTot2());
            cluster_tot2_from_eta2_layers[cluster.getLayer()].push_back(cluster.getTot2());
        }
    }
}

// Utility functions to push track-level data into the corresponding vectors for tree filling
void DataProcesser::pushBackTrackData(const Track& track) {
    if (track.getSide() == Track::ETA1) {
        track_length_eta1.push_back(track.getLayerCount());
        track_width_eta1.push_back(track.getWidth());
        track_size_eta1.push_back(track.getSize());
    } else if (track.getSide() == Track::ETA2) {
        track_length_eta2.push_back(track.getLayerCount());
        track_width_eta2.push_back(track.getWidth());
        track_size_eta2.push_back(track.getSize());
    }
}

// ==========================================================================================
// Main processing pipeline functions
// ==========================================================================================

// ------------------------------------------------------------------------------------------
// First stage: Entry point of the processing pipeline: Processing txt file and fill InputData tree with raw hit data
void DataProcesser::processInputData() {

    namespace fs = std::filesystem;

    if (!fs::exists(_input_path)) {
        std::cerr << "ERROR: Input path does not exist: " << _input_path << std::endl;
        return;
    }

    std::cout << "Processing file: " << _input_path << std::endl;
    if (_format == InputFormat::DecodedWords) {
        // For now do nothing since all the files are in the new filedump format
    } else {
        // Firts pass to fill InputData tree with raw hit information
        processDataFiledump(_input_path);

        // Second pass after InputData tree is filled and background rejection is applied
        processDataInputTree(_output_file);
    }
}

// Function to handle raw hit extraction (for filling InputData tree).
void DataProcesser::processDataFiledump(const std::string& file_path) {
    // Open the file and read line by line, extracting raw hit data in word format and filling the InputData tree
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open file: " << file_path << std::endl;
        return;
    }

    std::string token;
    bool in_event = false;
    int clk = 0;
    std::vector<int> event_clks;
    std::vector<int> event_words;

    while (infile >> token) {
        if (token.length() <= 512) {
            continue;
        }

        std::string header = token.substr(0, 2);
        if (header == "0a") {
            // PASS 2: Loop over all so far accumulated words from this event and fill the InputData tree with raw hit information
            if (in_event && !event_words.empty()) {

                // Extract raw hit information from the words and push back into raw data vectors
                for (size_t i = 0; i < event_words.size(); ++i) {
                    processSingleWord(event_clks[i], event_words[i]);
                }

                // After processing all words for this event, fill the InputData tree with the accumulated raw hit data
                if (event_words.size() > 0) {
                    _input_data_tree->Fill();    // Fill the tree with raw hit data for this event after processing all words
                }
            }

            // Prepare for the next event
            in_event = true;        // Set flag to indicate we are now processing an event
            clk = 0;                // Reset clock counter at the start of each event
            event_clks.clear();     // Clear event clocks for the new event
            event_words.clear();    // Clear event words for the new event
            clearRawDataVectors();  // Clear raw data vectors for the new event
        }

        if (!in_event) continue; // Skip any words until we encounter the first "0a" header indicating the start of an event

        // PASS 1: Extract words for this event and accumulate them in event_words and event_clks
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

            if (!valid_hex) continue;

            int dct = (word >> 28) & 0xF;
            if (dct == 1) clk++;

            if ((word & 0x0FFFFFFF) == 0x05555555) continue;

            event_clks.push_back(clk);
            event_words.push_back(static_cast<int>(word));
        }
    }

    // Process any remaining event at the end of file
    if (!event_words.empty()) {
        // PASS 2.5: Process the last event if the file does not end with a new event header
        for (size_t i = 0; i < event_words.size(); ++i) {
            processSingleWord(event_clks[i], event_words[i]);
        }
        if (event_words.size() > 0) {
            _input_data_tree->Fill();
        }
    }

    infile.close();
}

// Process a single word and accumulate its data into the current event's hit vector
void DataProcesser::processSingleWord(int clk, int word) {
    if (word == EMPTY_WORD) return;

    // Extract raw hit information from the word and fill the DCTWord struct
    _dct_word.clk = clk;
    decodeDCTWord(word);

    // Accumulate raw hit data into vector
    pushBackWordData(_dct_word);
}

// Utility function for extracting information from the raw DCT word
void DataProcesser::decodeDCTWord(int word) {
    _dct_word.channel = (word >> 20) & 0xFF;          // Bits 20-27
    _dct_word.rise = word & 0x01;                     // Bit 0

    if (_dct_word.rise == 1) {                      // Rising edge
        _dct_word.raw_bcid = word >> 12 & 0xFF;       // Bits 12-19
        _dct_word.raw_time_eta1 = word >> 7 & 0x1F;   // Bits 7-11
        _dct_word.raw_time_eta2 = word >> 1 & 0x1F;   // Bits 1-5

    } else {                                        // Falling edge
        _dct_word.raw_bcid = word >> 11 & 0x1FF;      // Bits 11-19
        _dct_word.raw_time_eta1 = word >> 6 & 0x1F;   // Bits 6-10
        _dct_word.raw_time_eta2 = word >> 1 & 0x1F;   // Bits 1-5
    }
}

// ------------------------------------------------------------------------------------------
// Run the main processing loop over the InputData tree, recreate Hit objects and process events
void DataProcesser::processDataInputTree(TFile* root_file) {

    // Check if the InputData tree is filled and can be accessed
    if (!root_file || !root_file->Get("InputData")) {
        std::cerr << "ERROR: Input file is null or InputData tree not found." << std::endl;
        return;
    }
    // Cast to TTree* to access tree-specific methods
    TTree* tree = dynamic_cast<TTree*>(root_file->Get("InputData"));
    if (!tree) {
        std::cerr << "ERROR: InputData is not a TTree." << std::endl;
        return;
    }
    int n_entries = tree->GetEntries();
    if (n_entries == 0) {
        std::cerr << "ERROR: InputData tree is empty." << std::endl;
        return;
    }

    // Read back raw data and recreate Hit objects
    for (int i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);

        // Accumulate hits for the current event
        for (size_t j = 0; j < hit_clk.size(); ++j) {

            // Recreate a Hit object from the raw hit data
            Hit hit(hit_clk[j], hit_channel[j], hit_raw_bcid[j], hit_raw_time1[j], hit_raw_time2[j], hit_rise[j]);
            hit.setIdx(current_event_hits.size());
            current_event_hits.push_back(hit);
        }

        if (current_event_hits.empty()) continue;   // Skip event if there are no hits after background rejection

        // After accumulating hits for this event, find BC0 (smallest raw BCID among hits in this event) and apply BCID correction to all hits in this event
        int min_raw_bcid = 256;
        for (const auto& hit : current_event_hits) {
            if (hit.getRawBCID() % 256 < min_raw_bcid) {
                min_raw_bcid = hit.getRawBCID() % 256;
            }
        }
        BC0 = min_raw_bcid;
        for (auto& hit : current_event_hits) {
            hit.applyBCIDCorrection(BC0);
        }

        // Now we can call processEvent to process this event
        processEvent(efficiency_counters, efficiency_counters_tracks);
        current_event_number++;
        current_event_hits.clear();
        clearEventVectors();
    }

    // Once all counters are updated, calculate the final efficiencies and create respective histograms
    updateEfficiencies();
    createHistograms();
}

/*
/// OLD:
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
        processSingleWord(clk, word, efficiency_counters, efficiency_counters_tracks);
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
*/

// Process a complete event that has been accumulated
void DataProcesser::processEvent(EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks) {
    n_hits = current_event_hits.size();
    // std::cout << "\n" << std::endl;
    // std::cout << "╔═════════════════════════════════════════════════════════════════╗" << std::endl;
    // std::cout << "║ Processing event No. " << current_event_number << " with " << n_hits << " hits" << std::endl;

    // Create an Event object and transfer hits using move semantics, passing the counters reference
    Event event(current_event_number, std::move(current_event_hits), counters, counters_tracks, _use_external_trigger);

    // Extract trigger information from the event
    event.extractTriggerTime();

    // Calculate ToT before clusterization
    event.calculateTOT();

    // Push back processed hit data into vectors for tree filling
    pushBackProcessedData(event);

    // Fill the processed data tree
    if (_processed_data_tree && proc_layer.size() > 0) {
        _processed_data_tree->Fill();
    }

    // Clusterization
    event.clusterize();

    // Update hit cluster IDs after clusterization
    updateClusterIDs(event);

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
    if (_clusterization_tree && (cluster_size_eta1.size() > 0 || cluster_size_eta2.size() > 0)) {
        _clusterization_tree->Fill();
    }

    // Track reconstruction
    event.reconstructTracks();

    // Update track IDs after track reconstruction
    updateTrackIDs(event);

    // Check if hits are in valid tracks
    isHitInValidTrack(event);

    // Push back track-level data for both eta1 and eta2 sides
    for (const auto& track : event.getTracks()) {
        pushBackTrackData(track);
    }

    // Reset and update efficiency flags for each new event
    event.updateEfficiencyFlags(_dt_max, _dt_min);

    // Update efficiency counters based on reconstructed tracks
    event.updateEfficiencyCounters();

    // Fill the track reconstruction tree (write if either eta1 or eta2 has tracks)
    if (_track_reconstruction_tree && (track_length_eta1.size() > 0 || track_length_eta2.size() > 0)) {
        _track_reconstruction_tree->Fill();
    }
}

// Utility function to update hit cluster IDs after clusterization
void DataProcesser::updateClusterIDs(const Event& event) {
    for (const auto& cluster : event.getClustersEta1()) {
        for (const auto& hit: cluster.getHits()) {
            _cluster_id_from_eta1.push_back(hit->getClusterIDEta1());
        }
    }
    for (const auto& cluster : event.getClustersEta2()) {
        for (const auto& hit: cluster.getHits()) {
            _cluster_id_from_eta2.push_back(hit->getClusterIDEta2());
        }
    }
}

// Utility function to update hit track IDs after track reconstruction
void DataProcesser::updateTrackIDs(const Event& event) {
    for (const auto& track : event.getTracksEta1()) {
        for (const auto& hit: track.getHits()) {
            _track_id_from_eta1.push_back(hit->getTrackIDEta1());
        }
    }
    for (const auto& track : event.getTracksEta2()) {
        for (const auto& hit: track.getHits()) {
            _track_id_from_eta2.push_back(hit->getTrackIDEta2());
        }
    }
}

// NOTE: Not really needed as this can be done directly with track ID from the hit together with isValidTrack()
// Utility function to check if a hit is a part of a valid track based on its track ID
void DataProcesser::isHitInValidTrack(const Event& event) {

    std::unordered_map<int, bool> valid_eta1, valid_eta2;
    for (const auto& track : event.getTracksEta1()) {
        valid_eta1[track.getTrackID()] = track.isValidTrack();
    }
    for (const auto& track : event.getTracksEta2()) {
        valid_eta2[track.getTrackID()] = track.isValidTrack();
    }
    
    // Look up each hit
    for (const auto& hit : event.getHits()) {
        int id1 = hit.getTrackIDEta1();
        _in_valid_track_eta1.push_back(id1 != -1 && valid_eta1[id1]);
        
        int id2 = hit.getTrackIDEta2();
        _in_valid_track_eta2.push_back(id2 != -1 && valid_eta2[id2]);
    }
}