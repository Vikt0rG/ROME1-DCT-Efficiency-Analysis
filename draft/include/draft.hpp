#pragma once

#include <vector>
#include <map>
#include <cmath>

// ============================================================
// Hit Class: Represents a single detector hit
// ============================================================
class Hit {
private:
    // Raw data from the DCT
    int clk, channel, raw_bcid, raw_bcout;
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
    Hit(int clk, int word, int raw_bcout);

    // Accessors
    int getClk() const { return clk; }
    int getChannel() const { return channel; }
    int getLayer() const { return layer; }
    int getStrip() const { return strip; }
    int getTimeEta1() const { return time_eta1; }
    int getTimeEta2() const { return time_eta2; }
    int getBCID() const { return bcid; }
    int getBCOut() const { return bcout; }
    int getRawTimeEta1() const { return raw_time_eta1; }
    int getRawTimeEta2() const { return raw_time_eta2; }
    int getRise() const { return rise; }

    // Processing methods
    void decodeDCTWord(int word);
    void applyBCWrapAround(int BC0);
    void geometryMapping();       // Maps channel to a layer and strip

    // Predicates for clustering and track reconstruction
    bool inCluster();

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

    // Clusterization and global parameters
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
// DataAnalysis Class: Main analysis section
// ============================================================
class DataAnalysis {
private:
    std::vector<Event> events;

    // Efficiency counters
    int triggered_events_external;
    int triggered_events_rpc[3];

    int eta1_efficiency_counter[3];
    int eta2_efficiency_counter[3];
    int eta_or_efficiency_counter[3];
    int eta_and_efficiency_counter[3];
    
public:
    // Constructor
    DataAnalysis();

    // Main analysis pipeline
    void processRawData(const std::string& input_file);
    void analyzeEfficiency();
    void writeResults(const std::string& output_file);

    // Accessors
    int getEventCount() const { return events.size(); }
    const std::vector<Event>& getEvents() const { return events; }
};

