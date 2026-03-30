#pragma once

#include <vector>
#include <map>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include <utils/utils.hpp>
#include <types.hpp>

// ============================================================
// Hit Class: Represents a single detector hit
// ============================================================
class Hit {
private:
    // Raw data from the DCT
    int clk, channel, raw_bcid, raw_bcout, bcout;
    int raw_time_eta1, raw_time_eta2, rise;

    // Processed data
    int bcid, bcout;              // After BC wrap-around correction
    int layer, strip;             // Geometry info
    int time_eta1, time_eta2;     // Converted times

    // IDs (WIP: Decide if neeeded at all)
    int hit_id;                   // Unique identifier for this hit
    int cluster_id;               // ID of cluster this hit belongs to (if any)
    int track_id;                 // ID of track this hit belongs to (if any)

public:
    // Constructor from raw word
    Hit(int clk, int word, int bcout, int BC0);

    // Accessors
    int getClk() const { return clk; }
    int getChannel() const { return channel; }
    int getLayer() const { return layer; }
    int getStrip() const { return strip; }
    int getTimeEta1() const { return time_eta1; }
    int getTimeEta2() const { return time_eta2; }
    int getBCID() const { return bcid; }
    int getBCOut() const { return bcout; }
    int getRawBCID() const { return raw_bcid; }
    int getRawBCOut() const { return raw_bcout; }
    int getRawTimeEta1() const { return raw_time_eta1; }
    int getRawTimeEta2() const { return raw_time_eta2; }
    int getRise() const { return rise; }

    // Processing methods
    void decodeDCTWord(int word);
    void applyBCWrapAround(int BC0);
    void geometryMapping();       // Maps channel to a layer and strip

    // Predicates for clustering and track reconstruction
    bool inCluster();
    bool inTrack();

    // Query methods
    bool hasEta1Time() const { return time_eta1 != -1; }
    bool hasEta2Time() const { return time_eta2 != -1; }
    bool hasTimeInfo() const { return time_eta1 != -1 || time_eta2 != -1; }
};


// ============================================================
// Cluster Class: Group of consecutive hits in space/time
// ============================================================
class Cluster {
private:
    int cluster_id;               // Unique identifier for the cluster
    std::vector<Hit*> hits;       // Hits in this cluster (NOT copied, referenced)
    int center_hit_index;         // Index of first-in-time hit (cluster center)
    int tot1, tot2;               // Time-over-threshold for each coordinate

    static constexpr int TIME_WINDOW_TICKS = 18;  // 15 ns clustering time window
    static constexpr int MAX_STRIP_DISTANCE = 1;  // Only consecutive strips

public:
    // Constructor
    Cluster(Hit* first_hit);

    // Core clustering logic
    bool addHit(Hit* hit, int time_window_ticks = TIME_WINDOW_TICKS);

    // Accessors
    Hit* getCenterHit() const { return hits[center_hit_index]; }
    int getSize() const { return hits.size(); }
    int getStrip() const { return getCenterHit()->getStrip(); }
    int getLayer() const { return getCenterHit()->getLayer(); }
    int getCenterTime() const { return getCenterHit()->getTimeEta1(); }
    int getTot1() const { return tot1; }
    int getTot2() const { return tot2; }

    const std::vector<Hit*>& getHits() const { return hits; }

    // ToT calculation (pair rising/falling edges)
    void calculateToT();
};


// ============================================================
// Track Class: Hits aligned across detector layers
// ============================================================
class Track {
private:
    std::vector<Hit*> hits;       // All hits in track (NOT copied, referenced)
    int track_id;

    // Layer occupancy flags (for efficiency calculations)
    bool eta1_layers[3];          // Layer 0, 1, 2 for η1 coordinate
    bool eta2_layers[3];          // Layer 0, 1, 2 for η2 coordinate

    // Track quality parameters
    static constexpr int MIN_TRACK_LENGTH = 2;    // Minimum hits required
    static constexpr int LAYER_COUNT = 3;

public:
    // Constructor
    Track(int id, Hit* first_hit);

    // Track reconstruction logic
    bool addHit(Hit* hit);        // Add hit if aligned with existing hits

    // Accessors
    int getId() const { return track_id; }
    int getLength() const { return hits.size(); }
    const std::vector<Hit*>& getHits() const { return hits; }

    // Layer occupancy queries (for efficiency calculations)
    bool hasEta1Layer(int layer) const;
    bool hasEta2Layer(int layer) const;
    bool hasAllEta1Layers() const;
    bool hasAllEta2Layers() const;
    bool hasAllLayers() const;    // Both eta1 and eta2 in all layers

    // Track quality methods
    bool isValidTrack() const { return hits.size() >= MIN_TRACK_LENGTH; }
    int getLayerCount(int coordinate) const;  // coordinate: 0=eta1, 1=eta2
};


// ============================================================
// Event Class: Container for hits, clusters, and tracks
// ============================================================
class Event {
private:
    int event_number;

    // Raw and processed hits
    std::vector<Hit> hits;

    // Clusters for each coordinate
    std::vector<Cluster> clusters_eta1;
    std::vector<Cluster> clusters_eta2;

    // Tracks for each coordinate
    std::vector<Track> tracks_eta1;
    std::vector<Track> tracks_eta2;

    // Trigger information
    int trigger_time;
    int trigger_channel;

    // Global parameters
    static constexpr int EMPTY_WORD = 0x5555555;
    static constexpr int CLUSTERING_TIME_WINDOW = 18;  // ticks
    static constexpr int TRIGGER_CHANNEL = 143;

public:
    // Constructor
    Event(int event_number);

    // Data input
    void addHit(int clk, int word, int raw_bcout);

    // Processing pipeline
    void clusterize();                    // Form clusters from hits
    void reconstructTracks();             // Form tracks from clusters/hits
    void calculateEfficiency();           // Efficiency counters and calculations

    // Extract trigger information
    void extractTriggerTime(int dt_min, int dt_max);

    // Accessors
    int getEventNumber() const { return event_number; }
    int getTriggerTime() const { return trigger_time; }
    int getHitCount() const { return hits.size(); }

    const std::vector<Hit>& getHits() const { return hits; }
    const std::vector<Cluster>& getClustersEta1() const { return clusters_eta1; }
    const std::vector<Cluster>& getClustersEta2() const { return clusters_eta2; }
    const std::vector<Track>& getTracksEta1() const { return tracks_eta1; }
    const std::vector<Track>& getTracksEta2() const { return tracks_eta2; }

    // Query methods for efficiency
    int getClusterCountEta1() const { return clusters_eta1.size(); }
    int getClusterCountEta2() const { return clusters_eta2.size(); }
    int getTrackCountEta1() const { return tracks_eta1.size(); }
    int getTrackCountEta2() const { return tracks_eta2.size(); }
};


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
    std::vector<int> track_length_eta1, track_length_eta2;

    // Event state management
    Event* current_event;
    int current_event_number;
    std::vector<Hit> current_event_hits;
    int BC0;  // BC0 reference for current event
    
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
    void processInputData(const std::string& input_path);
    void processFile(const std::string& file_path);
    void processSingleWord(int clk, int word, int raw_bcout);
    int extractRawBCID(int word);
    void pushBackHitData(const Hit& hit);
    void processEvent();
    
    // Analysis
    void calculateEfficiencies();
    void createHistograms();
    
    // Cleanup
    void clearEventVectors();
    void writeAndClose();
};

