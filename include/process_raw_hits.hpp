#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <vector>
#include <fstream>

// ============================================================================================================
// Data structure to hold efficiency counter data
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

// ============================================================================================================
// Setup functions

/// Create output ROOT file and define all trees with branches
void setupOutputFile(TFile*& output_file, 
                     TTree*& input_data_tree, 
                     TTree*& processed_data_tree, 
                     TTree*& clusterization_tree);

/// Setup all efficiency counter arrays to zero
void initializeEfficiencyCounters(EfficiencyCounters& counters);

// ============================================================================================================
// Word processing functions

/// Extract bit fields from DCT word (channel, rise, bcid, time1, time2)
void extractWordData(int word, int& channel, int& rise, int& raw_bcid, int& raw_time1, int& raw_time2);

/// Apply BC wrap-around correction to BCID values
void correctBCWrapAround(int raw_bcid, int raw_bcout, int BC0, int& bcid, int& bcout);

/// Store raw hit data in vectors
void storeRawHit(int clk, int word, int channel, int rise, int raw_bcid, int bcid, 
                 int raw_time1, int raw_time2, int raw_bcout, int bcout,
                 std::vector<int>& hit_clk, std::vector<int>& hit_channel, 
                 std::vector<int>& hit_raw_bcid, std::vector<int>& hit_bcid,
                 std::vector<int>& hit_time1, std::vector<int>& hit_time2, 
                 std::vector<int>& hit_rise, std::vector<int>& hit_raw_bcout, 
                 std::vector<int>& hit_bcout);

// ============================================================================================================
// Event processing functions

/// Extract trigger time from trigger hit (channel 143)
int extractTriggerTime(size_t n_hits, int trig_channel, 
                        const std::vector<int>& hit_channel, 
                        const std::vector<int>& hit_rise,
                        const std::vector<int>& hit_time2, 
                        const std::vector<int>& hit_bcid);

/// Process individual hit in event - efficiency counting and data storage
void processHitsInEvent(size_t idx_hit, int n_hits, int trigger_time,
                        int dt_max, int dt_min, int trig_channel,
                        const std::vector<int>& hit_channel,
                        const std::vector<int>& hit_rise,
                        const std::vector<int>& hit_time1,
                        const std::vector<int>& hit_time2,
                        const std::vector<int>& hit_bcid,
                        std::vector<int>& proc_layer,
                        std::vector<int>& proc_strip,
                        std::vector<int>& proc_time1,
                        std::vector<int>& proc_time2,
                        std::vector<int>& proc_dt_time1_time2,
                        std::vector<int>& proc_dt_time1_trigger,
                        std::vector<int>& proc_dt_time2_trigger,
                        bool* eta1_hit_flags,
                        bool* eta2_hit_flags);

/// Main event processing loop - called when clk==127
void processEvent(int& n_event, int n_hits,
                  std::vector<int>& hit_clk, std::vector<int>& hit_channel, 
                  std::vector<int>& hit_raw_bcid, std::vector<int>& hit_bcid,
                  std::vector<int>& hit_time1, std::vector<int>& hit_time2, 
                  std::vector<int>& hit_rise, std::vector<int>& hit_raw_bcout, 
                  std::vector<int>& hit_bcout,
                  std::vector<int>& proc_layer, std::vector<int>& proc_strip,
                  std::vector<int>& proc_time1, std::vector<int>& proc_time2,
                  std::vector<int>& proc_dt_time1_time2,
                  std::vector<int>& proc_trigger_time,
                  std::vector<int>& proc_dt_time1_trigger, std::vector<int>& proc_dt_time2_trigger,
                  std::vector<int>& proc_tot1, std::vector<int>& proc_tot2,
                  std::vector<int>& cluster_size1, std::vector<int>& cluster_size2,
                  std::vector<int>& cluster_tot1, std::vector<int>& cluster_tot2,
                  int dt_max, int dt_min, int trig_channel,
                  EfficiencyCounters& counters);

/// Clear all hit vectors for next event
void clearEventVectors(std::vector<int>& hit_clk, std::vector<int>& hit_channel,
                       std::vector<int>& hit_raw_bcid, std::vector<int>& hit_bcid,
                       std::vector<int>& hit_time1, std::vector<int>& hit_time2,
                       std::vector<int>& hit_rise, std::vector<int>& hit_raw_bcout,
                       std::vector<int>& hit_bcout, std::vector<int>& proc_layer,
                       std::vector<int>& proc_strip, std::vector<int>& proc_time1,
                       std::vector<int>& proc_time2, std::vector<int>& proc_dt_time1_time2,
                       std::vector<int>& proc_trigger_time,
                       std::vector<int>& proc_dt_time1_trigger, std::vector<int>& proc_dt_time2_trigger,
                       std::vector<int>& proc_tot1, std::vector<int>& proc_tot2,
                       std::vector<int>& cluster_size1, std::vector<int>& cluster_size2,
                       std::vector<int>& cluster_tot1, std::vector<int>& cluster_tot2);

// ============================================================================================================
// Efficiency calculation functions

/// Calculate efficiency values and errors from counters
void calculateEfficiencies(const EfficiencyCounters& counters,
                          double* eta1_eff_ext, double* eta2_eff_ext, 
                          double* eta_or_eff_ext, double* eta_and_eff_ext,
                          double* eta1_eff_ext_err, double* eta2_eff_ext_err, 
                          double* eta_or_eff_ext_err, double* eta_and_eff_ext_err,
                          double* eta1_eff_rpc, double* eta2_eff_rpc, 
                          double* eta_or_eff_rpc, double* eta_and_eff_rpc,
                          double* eta1_eff_rpc_err, double* eta2_eff_rpc_err, 
                          double* eta_or_eff_rpc_err, double* eta_and_eff_rpc_err);

/// Create and fill efficiency histograms in output file
void createEfficiencyHistograms(TFile* output_file,
                               const double* eta1_eff_ext, const double* eta2_eff_ext,
                               const double* eta_or_eff_ext, const double* eta_and_eff_ext,
                               const double* eta1_eff_ext_err, const double* eta2_eff_ext_err,
                               const double* eta_or_eff_ext_err, const double* eta_and_eff_ext_err,
                               const double* eta1_eff_rpc, const double* eta2_eff_rpc,
                               const double* eta_or_eff_rpc, const double* eta_and_eff_rpc,
                               const double* eta1_eff_rpc_err, const double* eta2_eff_rpc_err,
                               const double* eta_or_eff_rpc_err, const double* eta_and_eff_rpc_err);
