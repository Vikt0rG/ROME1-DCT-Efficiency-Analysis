#pragma once

#include <vector>


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