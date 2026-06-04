#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include "types.hpp"
#include "hit.hpp"
#include "event.hpp"
#include "cluster.hpp"
#include "track.hpp"
#include "constants.hpp"


// ============================================================
// DataProcesser Class: Main data processing section
// ============================================================
class DataProcesser {
private:
    TFile* _output_file = nullptr;

    // Trees
    TTree* input_data_tree = nullptr;
    TTree* processed_data_tree = nullptr;
    TTree* clusterization_tree = nullptr;
    TTree* track_reconstruction_tree = nullptr;

    // Efficiency counters and results
    EfficiencyCounters efficiency_counters{};
    EfficiencyCountersTracks efficiency_counters_tracks{};
    EfficiencyResults efficiency_results{};
    EfficiencyResultsTracks efficiency_results_tracks{};

    // Raw data vectors and structs
    DCTWord _dct_word;
    std::vector<int> hit_clk, hit_channel, hit_raw_bcid;
    std::vector<int> hit_raw_time1, hit_raw_time2, hit_rise;

    // Flag vector after background rejection
    std::vector<bool> is_signal;

    // Processed data vectors
    std::vector<int> proc_layer, proc_strip, proc_bc0, proc_bcid, proc_time1, proc_time2;
    std::vector<int> proc_dt_time1_time2, proc_trigger_time;
    std::vector<int> proc_dt_time1_trigger, proc_dt_time2_trigger;
    std::vector<int> proc_tot1, proc_tot2;

    // Cluster and track vectors
    std::vector<int> cluster_size_eta1, cluster_size_eta2;
    std::vector<int> cluster_tot1_from_eta1, cluster_tot2_from_eta1, cluster_tot1_from_eta2, cluster_tot2_from_eta2;    // Cluster center ToT from eta1 and eta2 clusterings
    std::array<std::vector<int>, 3> cluster_size_eta1_layers, cluster_size_eta2_layers;
    std::array<std::vector<int>, 3> cluster_tot1_from_eta1_layers, cluster_tot2_from_eta1_layers, cluster_tot1_from_eta2_layers, cluster_tot2_from_eta2_layers;  // Layer-specific cluster data
    std::vector<int> track_length_eta1, track_length_eta2, track_width_eta1, track_width_eta2, track_size_eta1, track_size_eta2;

    // ID vectors for hits in clusters and tracks
    std::vector<int> _cluster_id_from_eta1, _cluster_id_from_eta2, _track_id;

    // Event state management
    int BC0 = -100;  // BC0 reference for current event: Set to -100 as an invalid default value to detect if it was properly set
    int current_event_number = 0;
    int n_hits = 0;
    std::vector<Hit> current_event_hits;
    
    // Time window parameters for efficiency calculation
    int _dt_max = 20;
    int _dt_min = 0;

    // Flags
    bool _use_external_trigger = true;
    bool _reject_background = false;


public:
    enum class InputFormat {
        FiledumpPackets,
        DecodedWords
    };

    DataProcesser();
    ~DataProcesser();
    
    // Setup
    void setupOutputFile();
    void setupBranches();
    void initializeCounters();

    // Accessors for vectors
    std::vector<int>& getHitClk() { return hit_clk; }
    std::vector<int>& getHitChannel() { return hit_channel; }
    std::vector<int>& getHitRawBcid() { return hit_raw_bcid; }
    std::vector<int>& getHitRawTime1() { return hit_raw_time1; }
    std::vector<int>& getHitRawTime2() { return hit_raw_time2; }
    std::vector<int>& getHitRise() { return hit_rise; }

    std::vector<int>& getProcLayer() { return proc_layer; }
    std::vector<int>& getProcStrip() { return proc_strip; }
    std::vector<int>& getProcBC0() { return proc_bc0; }
    std::vector<int>& getProcBCID() { return proc_bcid; }
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
    std::vector<int>& getClusterTot1FromEta1() { return cluster_tot1_from_eta1; }
    std::vector<int>& getClusterTot2FromEta1() { return cluster_tot2_from_eta1; }
    std::vector<int>& getClusterTot1FromEta2() { return cluster_tot1_from_eta2; }
    std::vector<int>& getClusterTot2FromEta2() { return cluster_tot2_from_eta2; }

    std::vector<int>& getTrackLengthEta1() { return track_length_eta1; }
    std::vector<int>& getTrackLengthEta2() { return track_length_eta2; }
    std::vector<int>& getTrackWidthEta1() { return track_width_eta1; }
    std::vector<int>& getTrackWidthEta2() { return track_width_eta2; }
    std::vector<int>& getTrackSizeEta1() { return track_size_eta1; }
    std::vector<int>& getTrackSizeEta2() { return track_size_eta2; }

    // Processing
    void processInputData(const std::string& input_path, const int dt_max, const int dt_min, InputFormat format = InputFormat::FiledumpPackets, bool use_external_trigger_arg = true);

    // NEW PIPELINE
    void processInputData(const std::string& input_path, const int dt_max, const int dt_min, InputFormat format, bool use_external_trigger_arg, bool reject_background_arg);
    void processDataFiledump(const std::string& file_path);
    void processSingleWord(int clk, int word);
    void decodeDCTWord(int word);
    void applyBackgroundRejection();
    void processDataInputTree(TFile* root_file);
    void processEvent(EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks);
    void updateClusterIDs(const Event& event);

    void pushBackWordData(const DCTWord&);
    void pushBackProcessedData(const Event&);
    void pushBackClusterData(const Cluster&);
    void pushBackTrackDataEta1(const Track&);
    void pushBackTrackDataEta2(const Track&);

    // Analysis
    void updateEfficiencies();
    void createHistograms();

    // Cleanup
    void clearRawDataVectors();
    void clearEventVectors();
};