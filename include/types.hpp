#pragma once

#include <string>
#include <vector>
#include <array>


// ==========================================================================================
// Common types and structures used across the project

/// @struct DCTWord
/// @brief Represents a single DCT word with its decoded fields.
/// @param clk Clock cycle at which the word exits the DCT.
/// @param channel Channel number of the hit.
/// @param raw_bcid Raw BCID of the hit.
/// @param raw_time_eta1 Raw hit time for η1 side.
/// @param raw_time_eta2 Raw hit time for η2 side.
/// @param rise Hit's rising edge flag.
struct DCTWord {
    int clk = 0;
    int channel = 0;
    int raw_bcid = 0;
    int raw_time_eta1 = 0;
    int raw_time_eta2 = 0;
    int rise = 0;
};

/// @struct SelectionMask
/// @brief Represents a selection mask with a name and a boolean vector indicating selected entries.
/// @param name Name of the selection mask.
/// @param mask Boolean vector where true indicates selected entries.
struct SelectionMask {
    std::string name;
    std::vector<bool> mask;

    SelectionMask(const std::string& n, std::vector<bool>&& m) 
        : name(n), mask(std::move(m)) {}
};

// ==========================================================================================
// Efficiency calculation related structures
struct EfficiencyFlags {
    bool eta1_layer[3] = {false, false, false};         // Layer 0, 1, 2 for η1 side
    bool eta2_layer[3] = {false, false, false};         // Layer 0, 1, 2 for η2 side
    bool eta1_layer_track[3] = {false, false, false};   // Whether layer has a hit in a valid track for η1
    bool eta2_layer_track[3] = {false, false, false};   // Whether layer has a hit in a valid track for η2
};

struct EfficiencyCounters {
    int triggered_events_external;
    int triggered_events_rpc[3];
    
    int eta1_efficiency_counter[3];
    int eta2_efficiency_counter[3];
    int eta_or_efficiency_counter[3];
    int eta_and_efficiency_counter[3];
    
    int eta1_efficiency_counter_rpc[3];
    int eta2_efficiency_counter_rpc[3];
    int eta_or_efficiency_counter_rpc[3];
    int eta_and_efficiency_counter_rpc[3];
};

struct EfficiencyCountersTracks {
    int track_triggered_events_external;
    int track_triggered_events_rpc[3];
    
    int track_eta1_efficiency_counter[3];
    int track_eta2_efficiency_counter[3];
    int track_eta_or_efficiency_counter[3];
    int track_eta_and_efficiency_counter[3];
    
    int track_eta1_efficiency_counter_rpc[3];
    int track_eta2_efficiency_counter_rpc[3];
    int track_eta_or_efficiency_counter_rpc[3];
    int track_eta_and_efficiency_counter_rpc[3];
};

struct EfficiencyResults {
    double eta1_efficiency_external[3];
    double eta2_efficiency_external[3];
    double eta_or_efficiency_external[3];
    double eta_and_efficiency_external[3];
    
    double eta1_efficiency_external_error[3];
    double eta2_efficiency_external_error[3];
    double eta_or_efficiency_external_error[3];
    double eta_and_efficiency_external_error[3];

    double eta1_efficiency_rpc[3];
    double eta2_efficiency_rpc[3];
    double eta_or_efficiency_rpc[3];
    double eta_and_efficiency_rpc[3];

    double eta1_efficiency_rpc_error[3];
    double eta2_efficiency_rpc_error[3];
    double eta_or_efficiency_rpc_error[3];
    double eta_and_efficiency_rpc_error[3];
};

struct EfficiencyResultsTracks {
    double track_eta1_efficiency_external[3];
    double track_eta2_efficiency_external[3];
    double track_eta_or_efficiency_external[3];
    double track_eta_and_efficiency_external[3];

    double track_eta1_efficiency_external_error[3];
    double track_eta2_efficiency_external_error[3];
    double track_eta_or_efficiency_external_error[3];
    double track_eta_and_efficiency_external_error[3];

    double track_eta1_efficiency_rpc[3];
    double track_eta2_efficiency_rpc[3];
    double track_eta_or_efficiency_rpc[3];
    double track_eta_and_efficiency_rpc[3];

    double track_eta1_efficiency_rpc_error[3];
    double track_eta2_efficiency_rpc_error[3];
    double track_eta_or_efficiency_rpc_error[3];
    double track_eta_and_efficiency_rpc_error[3];
};

// ==========================================================================================
// Measurement configuration and summary structures for the measurement analysis part
struct MeasurementEntry {
    std::string name;
    std::string measurement_type;
    std::string root_file;
    std::string mixture;
    std::string source;
    double filter = 0.0;
    int lv_setting = 0;
    double hv = 0.0;
    double other_hv = 0.0;
    double scanned_hv = 0.0;
    int layer = 0;
};

struct PerFileStats {
    std::string name;
    int layer = 0;
    double hv = 0.0;

    EfficiencyResults efficiency_results;
    EfficiencyResultsTracks efficiency_results_tracks;

    double avg_cluster_size_eta1 = 0.0;
    double avg_cluster_size_eta2 = 0.0;
    double avg_cluster_size_eta1_layers[3] = {0.0, 0.0, 0.0};
    double avg_cluster_size_eta2_layers[3] = {0.0, 0.0, 0.0};
    double noise_rate = 0.0;
    double noise_rate_eta1[3] = {0.0, 0.0, 0.0};
    double noise_rate_eta2[3] = {0.0, 0.0, 0.0};
};

struct SummaryStats {
    std::string config_path;
    std::string measurement_type;
    std::vector<MeasurementEntry> entries;
    std::vector<PerFileStats> per_file_stats;
};