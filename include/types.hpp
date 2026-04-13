#pragma once

#include <vector>


struct EfficiencyFlags {
    bool eta1_layer[3];          // Layer 0, 1, 2 for η1 side
    bool eta2_layer[3];          // Layer 0, 1, 2 for η2 side
    bool eta1_layer_track[3];    // Whether layer has a hit in a valid track for η1
    bool eta2_layer_track[3];    // Whether layer has a hit in a valid track for η2
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