#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "types.hpp"
#include "hit.hpp"
#include "event.hpp"
#include "cluster.hpp"
#include "track.hpp"


// ============================================================
// DataAnalysis Class: Main analysis section: WIP
// ============================================================
class DataAnalyzer {
private:
    TFile* output_file;
    TTree* input_data_tree;
    TTree* processed_data_tree;
    TTree* clusterization_tree;
    TTree* track_reconstruction_tree;
    EfficiencyCounters efficiency_counters;
    EfficiencyCountersTracks efficiency_counters_tracks;
    EfficiencyResults efficiency_results;
    EfficiencyResultsTracks efficiency_results_tracks;
    
    // Raw data vectors
    std::vector<int> hit_clk, hit_channel, hit_raw_bcid, hit_bcid;
    std::vector<int> hit_time1, hit_time2, hit_rise, hit_raw_bcout, hit_bcout;
    
    // Processed data vectors
    std::vector<int> proc_layer, proc_strip, proc_time1, proc_time2;
    std::vector<int> proc_dt_time1_time2, proc_trigger_time;
    std::vector<int> proc_dt_time1_trigger, proc_dt_time2_trigger;
    std::vector<int> proc_tot1, proc_tot2;
    
    // Cluster and track vectors
    std::vector<int> cluster_size_eta1, cluster_size_eta2, cluster_tot1, cluster_tot2;
    std::vector<int> track_length_eta1, track_length_eta2, track_width_eta1, track_width_eta2, track_size_eta1, track_size_eta2;

    // Event state management
    Event* current_event;
    int current_event_number;
    int n_hits;
    std::vector<Hit> current_event_hits;
    int BC0;  // BC0 reference for current event
    
    // Time window parameters for efficiency calculation
    int dt_max;
    int dt_min;
    
    // Constants
    static constexpr int EMPTY_WORD = 0x5555555;  // Pattern indicating empty/skipped word

public:
    DataAnalyzer();
    ~DataAnalyzer();
    
    // Setup
    void setupOutputFile();
    void setupBranches();
    void initializeCounters();
    
    // Accessors for vectors
    std::vector<int>& getHitClk() { return hit_clk; }
    std::vector<int>& getHitChannel() { return hit_channel; }
    std::vector<int>& getHitRawBcid() { return hit_raw_bcid; }
    std::vector<int>& getHitBcid() { return hit_bcid; }
    std::vector<int>& getHitTime1() { return hit_time1; }
    std::vector<int>& getHitTime2() { return hit_time2; }
    std::vector<int>& getHitRise() { return hit_rise; }
    std::vector<int>& getHitRawBcout() { return hit_raw_bcout; }
    std::vector<int>& getHitBcout() { return hit_bcout; }
    
    std::vector<int>& getProcLayer() { return proc_layer; }
    std::vector<int>& getProcStrip() { return proc_strip; }
    std::vector<int>& getProcTime1() { return proc_time1; }
    std::vector<int>& getProcTime2() { return proc_time2; }
    std::vector<int>& getProcDtTime1Time2() { return proc_dt_time1_time2; }
    std::vector<int>& getProcTriggerTime() { return proc_trigger_time; }
    std::vector<int>& getProcDtTime1Trigger() { return proc_dt_time1_trigger; }
    std::vector<int>& getProcDtTime2Trigger() { return proc_dt_time2_trigger; }
    std::vector<int>& getProcTot1() { return proc_tot1; }
    std::vector<int>& getProcTot2() { return proc_tot2; }
    
    std::vector<int>& getClusterSizeEta1() { return cluster_size_eta1; }
    std::vector<int>& getClusterSizeEta2() { return cluster_size_eta2; }
    std::vector<int>& getClusterTot1() { return cluster_tot1; }
    std::vector<int>& getClusterTot2() { return cluster_tot2; }
    
    std::vector<int>& getTrackLengthEta1() { return track_length_eta1; }
    std::vector<int>& getTrackLengthEta2() { return track_length_eta2; }
    
    // Processing
    void processInputData(const std::string& input_path, const int dt_max, const int dt_min);
    void processFile(const std::string& file_path);
    void processSingleWord(int clk, int word, int raw_bcout, EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks);
    int extractRawBCID(int word);
    void pushBackHitData(const Hit& hit);
    void pushBackProcessedData(const Event& event);
    void pushBackClusterDataEta1(const Cluster& cluster);
    void pushBackClusterDataEta2(const Cluster& cluster);
    void pushBackTrackDataEta1(const Track& track);
    void pushBackTrackDataEta2(const Track& track);
    void processEvent(EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks);
    
    // Analysis
    void updateEfficiencies();
    void createHistograms();
    
    // Cleanup
    void clearEventVectors();
    void writeAndClose();
};