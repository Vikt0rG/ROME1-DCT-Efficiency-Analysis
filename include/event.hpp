#pragma once

#include <iostream>
#include "types.hpp"
#include "cluster.hpp"
#include "track.hpp"
#include "constants.hpp"


// ============================================================
// Event Class: Container for hits, clusters, and tracks
// ============================================================
class Event {
private:
    // Event metadata
    int event_number;

    // Trigger information
    int trigger_time;
    int trigger_channel;

    // Raw and processed hits
    std::vector<Hit> hits;

    // Cluster vectors
    std::vector<Cluster> clusters_eta1;
    std::vector<Cluster> clusters_eta2;

    // Track vectors
    std::vector<Track> tracks_eta1;
    std::vector<Track> tracks_eta2;

    // Efficiency flags and counters
    EfficiencyFlags efficiency_flags;

    // Reference to shared efficiency counters (for updating throughout event lifetime)
    EfficiencyCounters& efficiency_counters;
    EfficiencyCountersTracks& efficiency_counters_tracks;

    bool use_external_trigger;

public:
    // Constructor
    Event(int event_number, EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks, bool use_external_trigger);
    Event(int event_number, std::vector<Hit>&& hits_in, EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks, bool use_external_trigger);  // Move constructor for hits
    ~Event();  // Destructor

    // Extract trigger information
    void extractTriggerTime();

    // Processing pipeline
    void calculateTOT();                  // Calculate Time-over-Threshold for hits
    void clusterize();                    // Form clusters from hits
    void calculateTOTCluster();           // Calculate Time-over-Threshold for clusters
    void reconstructTracks();             // Form tracks from clusters/hits
    void updateEfficiencyFlags(
        const int dt_max,
        const int dt_min);                // Update efficiency flags based on time window
    void updateEfficiencyCounters();      // Update efficiency counters based on reconstructed tracks

    // Accessors
    int getEventNumber() const { return event_number; }
    int getTriggerChannel() const { return trigger_channel; }
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