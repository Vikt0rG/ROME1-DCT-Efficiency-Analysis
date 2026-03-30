#pragma once

#include "types.hpp"
#include "cluster.hpp"
#include "track.hpp"


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