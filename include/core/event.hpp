#pragma once

#include "core/types.hpp"
#include "core/cluster.hpp"
#include "core/track.hpp"
#include "core/constants.hpp"

// ============================================================
// Event Class: Container for hits, clusters, and tracks
// ============================================================
class Event {
private:
    // Event metadata
    int _event_number;

    // Trigger information
    int _trigger_time = -1;  // Trigger time initialized to invalid
    int _trigger_channel = TRIGGER_CHANNEL;  // Default trigger channel from constants

    // Raw and processed hits
    std::vector<Hit> _hits;

    // Cluster vectors
    std::vector<Cluster> _clusters;
    std::vector<Cluster> _clusters_eta1;
    std::vector<Cluster> _clusters_eta2;

    // Track vectors
    std::vector<Track> _tracks;
    std::vector<Track> _tracks_eta1;
    std::vector<Track> _tracks_eta2;

    // Efficiency flags and counters
    EfficiencyFlags _efficiency_flags;  // Automatically initialized to false for all layers and track flags

    // Reference to shared efficiency counters (for updating throughout event lifetime)
    EfficiencyCounters& _efficiency_counters;
    EfficiencyCounters& _efficiency_counters_tracks;

    bool _use_external_trigger;

public:
    // Constructor and destructor
    Event(int event_number, std::vector<Hit>&& hits_in, EfficiencyCounters& counters, EfficiencyCounters& counters_tracks, bool use_external_trigger);
    ~Event();

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

    // Utility methods
    std::pair<int, int> countEdges() const;  // Count edges in the event for efficiency calculations

    // Accessors
    int getEventNumber() const { return _event_number; }
    int getTriggerChannel() const { return _trigger_channel; }
    int getTriggerTime() const { return _trigger_time; }
    int getHitCount() const { return _hits.size(); }

    const std::vector<Hit>& getHits() const { return _hits; }
    const std::vector<Cluster>& getClusters() const { return _clusters; }
    const std::vector<Cluster>& getClustersEta1() const { return _clusters_eta1; }
    const std::vector<Cluster>& getClustersEta2() const { return _clusters_eta2; }
    const std::vector<Track>& getTracks() const { return _tracks; }
    const std::vector<Track>& getTracksEta1() const { return _tracks_eta1; }
    const std::vector<Track>& getTracksEta2() const { return _tracks_eta2; }

    // Query methods for efficiency
    int getClusterCountEta1() const { return _clusters_eta1.size(); }
    int getClusterCountEta2() const { return _clusters_eta2.size(); }
    int getTrackCountEta1() const { return _tracks_eta1.size(); }
    int getTrackCountEta2() const { return _tracks_eta2.size(); }
};