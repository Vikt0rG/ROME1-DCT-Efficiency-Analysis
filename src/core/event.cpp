#include <iostream>

#include "core/event.hpp"

// ============================================================
// Event class implementation
// ============================================================
// Constructor and destructor for Event class
Event::Event(int event_number, std::vector<Hit>&& hits_in, EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks, bool use_external_trigger) 
    : _event_number(event_number)
    , _hits(std::move(hits_in))
    , _efficiency_counters(counters)
    , _efficiency_counters_tracks(counters_tracks)
    , _use_external_trigger(use_external_trigger) {
}

// Destructor for the event class
Event::~Event() {
    // Cleanup code if necessary (currently not needed)
}

// Utility function to add a hit to extract trigger time
void Event::extractTriggerTime() {
    if (!_use_external_trigger) {
        _trigger_time = -1;
        return;
    }
    for (const auto& hit : _hits) {
        if (hit.getChannel() == _trigger_channel && hit.getRise() == 1) { // Check for rising edge on trigger channel
            _trigger_time = hit.getTimeEta2(); // Trigger signal is on the η2 side
            break; // Assuming only one trigger hit per event, we can break after finding it
        }
    }
}

// Utility function to calculate the number of rising and falling edges in the event
std::pair<int, int> Event::countEdges() const {
    int rising_count = 0;
    int falling_count = 0;
    for (const auto& hit : _hits) {
        if (hit.getRise() == 1) rising_count++;
        else if (hit.getRise() == 0) falling_count++;
    }
    return {rising_count, falling_count};
}

// Time Over Threshold calculation before clusterization
void Event::calculateTOT() {
    // Vectors to keep track which hits have been paired for ToT calculation to avoid double counting
    // Note: For pairing rely on eta1 side only unless it doesn't have time information, then use eta2 side.
    std::vector<bool> paired(_hits.size(), false);

    // Determine if we have more rising edges than falling edges in the event, which can indicate potential data issues (?)
    auto [rising_count, falling_count] = countEdges();
    bool more_rising_edges = rising_count > falling_count;
    bool more_falling_edges = !more_rising_edges;

    for (size_t i = 0; i < _hits.size(); i++) {
        Hit& hit = _hits[i];
        int tot1 = -1;
        int tot2 = -1;

        if (hit.getChannel() == _trigger_channel) continue; // Skip hits on trigger channel

        if (paired[i]) continue; // Skip already paired hits

        // Greedy approach: Skip hits with edges that are more common, to match those to the less common ones
        if ((hit.getRise() == 0 && more_falling_edges) || (hit.getRise() == 1 && more_rising_edges)) continue;

        // Find the closest in-time partner edge on the same channel with a different edge type
        int best_j = -1;
        int min_dt = INT_MAX;
        const int target_rise = (hit.getRise() == 1) ? 0 : 1;

        for (size_t j = 0; j < _hits.size(); j++) {
            Hit& potential_partner = _hits[j];
            if (&potential_partner == &hit || paired[j]) continue;              // Skip the same hit or already paired hits
            if (potential_partner.getChannel() != hit.getChannel()) continue;   // Skip hits from different channels
            if (potential_partner.getRise() != target_rise) continue;           // Skip hits with the same edge type

            // Use eta1 for matching if available, otherwise fallback to eta2
            int dt = -1;
            if (hit.hasEta1Time() && potential_partner.hasEta1Time()) {
                dt = potential_partner.getTimeEta1() - hit.getTimeEta1();
            } else if (hit.hasEta2Time() && potential_partner.hasEta2Time()) {
                dt = potential_partner.getTimeEta2() - hit.getTimeEta2();
            } else {
                continue; // Skip if neither hit has valid time information
            }

            if (dt > 0 && dt < min_dt) {
                min_dt = dt;
                best_j = static_cast<int>(j);
            }
        }

        // Now calculate both ToTs using the same best pair
        if (best_j != -1) {
            Hit& partner = _hits[best_j];
            paired[i] = true;
            paired[best_j] = true;

            // Calculate tot1 if both have eta1 times
            if (hit.hasEta1Time() && partner.hasEta1Time()) {
                tot1 = partner.getTimeEta1() - hit.getTimeEta1();
                if (tot1 > 0) {
                    hit.setTot1(tot1);
                    partner.setTot1(tot1);
                }
            }

            // Calculate tot2 if both have eta2 times
            if (hit.hasEta2Time() && partner.hasEta2Time()) {
                tot2 = partner.getTimeEta2() - hit.getTimeEta2();
                if (tot2 > 0) {
                    hit.setTot2(tot2);
                    partner.setTot2(tot2);
                }
            }
        }
    }
}

// Clusterization method to form clusters from hits
void Event::clusterize() {
    // Side η1 clusterization
    // Set cluster ID counter
    int cluster_id_counter_eta1 = 0;

    // Loop over the hits and form clusters based on time and strip adjacency
    for (Hit& hit : _hits) {

        // Initialize a new cluster object if the rising non-trigger hit is not already in a cluster and has valid time information for η1
        if (hit.getRise() == 0 || hit.getChannel() == _trigger_channel || hit.inClusterEta1() || hit.getTimeEta1() == -1) continue;
        // std::cout << "------------------------------------------------------------------" << std::endl;
        // std::cout << "Processing side η1: Hit: " << hit.getIdx() << "; Layer " << hit.getLayer() << "; Strip: " << hit.getStrip() << "; Time " << hit.getTimeEta1() << std::endl;
        Cluster cluster(&hit, Cluster::ETA1, hit.getLayer());
        cluster.setClusterID(cluster_id_counter_eta1);
        hit.setClusterIDEta1(cluster_id_counter_eta1);

        // Loop over the remaining hits to find cluster partners for the current hit
        for (Hit& potential_partner : _hits) {
            // Skip the same hit
            if (&potential_partner == &hit) continue;

            // Skip falling hits
            if (potential_partner.getRise() == 0) continue;

            // Skip hits from another layer
            if (potential_partner.getLayer() != hit.getLayer()) continue;

            // Skip hits that are already in a cluster
            if (potential_partner.inClusterEta1()) continue;

            // Skip hits without valid time information
            if (potential_partner.getTimeEta1() == -1) continue;

            // Check logic for cluster membership (time window and strip adjacency)
            if (cluster.addHit(&potential_partner)) potential_partner.setClusterIDEta1(cluster_id_counter_eta1);
        }

        // If no more cluster partners are found for the current hit, store cluster information and increment cluster ID counter for the next cluster
        _clusters.push_back(cluster);
        _clusters_eta1.push_back(cluster);
        cluster_id_counter_eta1++;

        /*
        // DEBUG: Print out parameters of the hits belonging to the cluster for the current hit
        std::cout << "******************************************************************" << std::endl;
        std::cout << "Resulting cluster: " << std::endl;
        for (Hit* cluster_hit : cluster.getHits()) {
            std::cout << cluster_hit->getIdx() << ": layer = " << cluster_hit->getLayer() << ", strip = " << cluster_hit->getStrip() << ", time = " << cluster_hit->getTimeEta1() << std::endl;
        }
        std::cout << "\n";
        */
    }

    // Side η2 clusterization (same logic as for η1)
    int cluster_id_counter_eta2 = 0;
    for (Hit& hit : _hits) {
        if (hit.getRise() == 0 || hit.getChannel() == _trigger_channel || hit.getTimeEta2() == -1 || hit.inClusterEta2()) continue;
        // std::cout << "------------------------------------------------------------------" << std::endl;
        // std::cout << "Processing side η2: Hit: " << hit.getIdx() << "; Layer " << hit.getLayer() << "; Strip: " << hit.getStrip() << "; Time " << hit.getTimeEta2() << std::endl;
        Cluster cluster(&hit, Cluster::ETA2, hit.getLayer());
        cluster.setClusterID(cluster_id_counter_eta2);
        hit.setClusterIDEta2(cluster_id_counter_eta2);

        for (Hit& potential_partner : _hits) {
            if (&potential_partner == &hit) continue;

            if (potential_partner.getRise() == 0) continue;

            if (potential_partner.getLayer() != hit.getLayer()) continue;

            if (potential_partner.inClusterEta2()) continue;

            if (potential_partner.getTimeEta2() == -1) continue;

            if (cluster.addHit(&potential_partner)) potential_partner.setClusterIDEta2(cluster_id_counter_eta2);
        }
        _clusters.push_back(cluster);
        _clusters_eta2.push_back(cluster);
        cluster_id_counter_eta2++;

        /*
        std::cout << "******************************************************************" << std::endl;
        std::cout << "Resulting cluster: " << std::endl;
        for (Hit* cluster_hit : cluster.getHits()) {
            std::cout << cluster_hit->getIdx() << ": layer = " << cluster_hit->getLayer() << ", strip = " << cluster_hit->getStrip() << ", time = " << cluster_hit->getTimeEta2() << std::endl;
        }
        std::cout << "\n";
        */
    }
}

// Time-over-threshold calculation for cluster centers (after clusterization)
void Event::calculateTOTCluster() {

    // Vectors to keep track which hits have been paired for ToT calculation to avoid double counting
    // Note: For pairing rely on eta1 side unless it doesn't have time information, then use eta2 side.
    std::vector<bool> paired_eta1(_hits.size(), false);
    std::vector<bool> paired_eta2(_hits.size(), false);

    // Check if there are more cluster centers for eta1/eta2 than falling edges
    auto [_, falling_count] = countEdges();
    bool more_eta1_centers = _clusters_eta1.size() > static_cast<size_t>(falling_count);
    bool more_eta2_centers = _clusters_eta2.size() > static_cast<size_t>(falling_count);

    /// Calculate ToT for cluster centers from eta1 clustering first
    for (size_t i = 0; i < _hits.size(); i++) {
        Hit& hit = _hits[i];

        if (hit.getChannel() == _trigger_channel) continue; // Skip hits on trigger channel

        if (paired_eta1[i]) continue; // Skip already paired hits

        // Greedy approach: Skip hits with edges that are more common, to match those to the less common ones
        if ((hit.getRise() == 0 && !more_eta1_centers) ||   // Edge is falling and those are more common, so skip it
            (hit.getRise() == 1 && more_eta1_centers)       // Edge is rising and those are more common, so skip it
        ) continue;

        // Find the closest in-time partner edge on the same channel with a different edge type for eta1 ToT calculation
        int best_j_eta1 = -1;
        int min_dt_eta1 = INT_MAX;
        const int target_rise = (hit.getRise() == 1) ? 0 : 1;

        for (size_t j = 0; j < _hits.size(); j++) {
            Hit& potential_partner = _hits[j];
            if (&potential_partner == &hit || paired_eta1[j]) continue;         // Skip the same hit or already paired hits
            if (potential_partner.getChannel() != hit.getChannel()) continue;   // Skip hits from different channels
            if (potential_partner.getRise() != target_rise) continue;           // Skip hits with the same edge type
            if (!hit.isClusterCenterEta1() && 
                !potential_partner.isClusterCenterEta1()
            ) continue; // Only consider pairs where at least one hit is a cluster center on eta1 side

            // Use eta1 for matching if available, otherwise fallback to eta2
            int dt = -1;
            if (hit.hasEta1Time() && potential_partner.hasEta1Time()) {
                dt = potential_partner.getTimeEta1() - hit.getTimeEta1();
            } else if (hit.hasEta2Time() && potential_partner.hasEta2Time()) {
                dt = potential_partner.getTimeEta2() - hit.getTimeEta2();
            } else {
                continue; // Skip if neither hit has valid time information
            }
            if (dt > 0 && dt < min_dt_eta1) {
                min_dt_eta1 = dt;
                best_j_eta1 = static_cast<int>(j);
            }
        }

        // Now calculate ToT for eta1 clustering if a partner is found
        if (best_j_eta1 != -1) {
            Hit& partner = _hits[best_j_eta1];
            paired_eta1[i] = true;
            paired_eta1[best_j_eta1] = true;

            // Calculate ToT1 if both have eta1 times
            if (hit.hasEta1Time() && partner.hasEta1Time()) {
                int tot1_from_eta1 = partner.getTimeEta1() - hit.getTimeEta1();
                // Get cluster ID from one of two hits
                int cluster_id = hit.isClusterCenterEta1() ? hit.getClusterIDEta1() : partner.getClusterIDEta1();
                Cluster& cluster = _clusters_eta1[cluster_id];
                cluster.setTot1(tot1_from_eta1);
            }

            // Calculate ToT2 if both have eta2 times
            if (hit.hasEta2Time() && partner.hasEta2Time()) {
                int tot2_from_eta1 = partner.getTimeEta2() - hit.getTimeEta2();
                // Get cluster ID from one of two hits
                int cluster_id = hit.isClusterCenterEta1() ? hit.getClusterIDEta1() : partner.getClusterIDEta1();
                Cluster& cluster = _clusters_eta1[cluster_id];
                cluster.setTot2(tot2_from_eta1);
            }
        }
    }

    // Second pass to calculate ToT for cluster centers from eta2 clustering
    for (size_t i = 0; i < _hits.size(); i++) {
        Hit& hit = _hits[i];

        if (hit.getChannel() == _trigger_channel) continue; // Skip hits on trigger channel

        if (paired_eta2[i]) continue; // Skip already paired hits

        // Greedy approach: Skip hits with edges that are more common, to match those to the less common ones
        if ((hit.getRise() == 0 && !more_eta2_centers) ||   // Edge is falling and those are more common, so skip it
            (hit.getRise() == 1 && more_eta2_centers)       // Edge is rising and those are more common, so skip it
        ) continue;

        // Find the closest in-time partner edge on the same channel with a different edge type for eta2 ToT calculation
        int best_j_eta2 = -1;
        int min_dt_eta2 = INT_MAX;
        const int target_rise = (hit.getRise() == 1) ? 0 : 1;

        for (size_t j = 0; j < _hits.size(); j++) {
            Hit& potential_partner = _hits[j];
            if (&potential_partner == &hit || paired_eta2[j]) continue;         // Skip the same hit or already paired hits
            if (potential_partner.getChannel() != hit.getChannel()) continue;   // Skip hits from different channels
            if (potential_partner.getRise() != target_rise) continue;           // Skip hits with the same edge type
            if (!hit.isClusterCenterEta2() && 
                !potential_partner.isClusterCenterEta2()
            ) continue; // Only consider pairs where at least one hit is a cluster center on eta2 side

            int dt = -1;
            if (hit.hasEta2Time() && potential_partner.hasEta2Time()) {
                dt = potential_partner.getTimeEta2() - hit.getTimeEta2();
            } else if (hit.hasEta1Time() && potential_partner.hasEta1Time()) {
                dt = potential_partner.getTimeEta1() - hit.getTimeEta1();
            } else {
                continue; // Skip if neither hit has valid time information
            }
            if (dt > 0 && dt < min_dt_eta2) {
                min_dt_eta2 = dt;
                best_j_eta2 = static_cast<int>(j);
            }
        }

        if (best_j_eta2 != -1) {
            Hit& partner = _hits[best_j_eta2];
            paired_eta2[i] = true;
            paired_eta2[best_j_eta2] = true;

            // Calculate ToT1 if both have eta1 times
            if (hit.hasEta1Time() && partner.hasEta1Time()) {
                int tot1_from_eta2 = partner.getTimeEta1() - hit.getTimeEta1();
                // Get cluster ID from one of two hits
                int cluster_id = hit.isClusterCenterEta2() ? hit.getClusterIDEta2() : partner.getClusterIDEta2();
                Cluster& cluster = _clusters_eta2[cluster_id];
                cluster.setTot1(tot1_from_eta2);
            }

            // Calculate ToT2 if both have eta2 times
            if (hit.hasEta2Time() && partner.hasEta2Time()) {
                int tot2_from_eta2 = partner.getTimeEta2() - hit.getTimeEta2();
                // Get cluster ID from one of two hits
                int cluster_id = hit.isClusterCenterEta2() ? hit.getClusterIDEta2() : partner.getClusterIDEta2();
                Cluster& cluster = _clusters_eta2[cluster_id];
                cluster.setTot2(tot2_from_eta2);
            }
        }
    }
}

// Track reconstruction from clusters
void Event::reconstructTracks() {
    // Side η1 track reconstruction logic
    int track_id_counter_eta1 = 0;

    // Loop over η1 cluster centers and form tracks based on time and strip alignment across layers
    for (Cluster& cluster : _clusters_eta1) {

        // Initialize a new track object if the cluster center is not already in a track
        // NOTE: No checks for rising edge, non-trigger channel, valid time information as these are already ensured for cluster centers at the clusterization level
        if (cluster.getCenterHit()->inTrackEta1()) continue;    // Skip if the cluster center is already in a track on eta1 side
        // std::cout << "------------------------------------------------------------------" << std::endl;
        // std::cout << "Processing track reconstruction for side η1: " << std::endl;
        Track track(cluster.getCenterHit(), Track::ETA1);
        track.setTrackID(track_id_counter_eta1);
        cluster.getCenterHit()->setTrackIDEta1(track_id_counter_eta1);

        // Loop over the remaining cluster centers to find track partners for the current cluster center
        for (Cluster& potential_partner_cluster : _clusters_eta1) {
            // Skip the same cluster
            if (&potential_partner_cluster == &cluster) continue;

            // Skip cluster centers from the same layer as the first cluster center
            if (potential_partner_cluster.getCenterHit()->getLayer() == cluster.getCenterHit()->getLayer()) continue;

            // Skip cluster centers that are already in a track
            if (potential_partner_cluster.getCenterHit()->inTrackEta1()) continue;

            // Check logic for track membership (strip window and time alignment)
            if (track.addHit(potential_partner_cluster.getCenterHit())) potential_partner_cluster.getCenterHit()->setTrackIDEta1(track_id_counter_eta1);
        }

        // If no more track partners are found for the current cluster center, store track information and increment track ID counter for the next track
        _tracks.push_back(track);
        _tracks_eta1.push_back(track);
        track_id_counter_eta1++;

        /*
        // DEBUG: Print out parameters of the hits belonging to the track for the current cluster center
        std::cout << "******************************************************************" << std::endl;
        std::cout << "Resulting track: " << std::endl;
        for (Hit* track_hit : track.getHits()) {
            std::cout << track_hit->getIdx() << ": layer = " << track_hit->getLayer() << ", strip = " << track_hit->getStrip() << ", time = " << track_hit->getTimeEta1() << std::endl;
            if (track.getLayerCount(0) == 0) std::exit(1); // Sanity check to ensure layer count is correctly calculated for eta1 tracks
        }
        std::cout << "\n";
        */
    }

    // Side η2 track reconstruction logic (same logic as for η1)
    int track_id_counter_eta2 = 0;
    for (Cluster& cluster : _clusters_eta2) {
        if (cluster.getCenterHit()->inTrackEta2()) continue;
        // std::cout << "------------------------------------------------------------------" << std::endl;
        // std::cout << "Processing track reconstruction for side η2: " << std::endl;
        Track track(cluster.getCenterHit(), Track::ETA2);
        track.setTrackID(track_id_counter_eta2);
        cluster.getCenterHit()->setTrackIDEta2(track_id_counter_eta2);

        for (Cluster& potential_partner_cluster : _clusters_eta2) {
            if (&potential_partner_cluster == &cluster) continue;

            if (potential_partner_cluster.getCenterHit()->getLayer() == cluster.getCenterHit()->getLayer()) continue;

            if (potential_partner_cluster.getCenterHit()->inTrackEta2()) continue;

            if (track.addHit(potential_partner_cluster.getCenterHit())) potential_partner_cluster.getCenterHit()->setTrackIDEta2(track_id_counter_eta2);
        }
        _tracks.push_back(track);
        _tracks_eta2.push_back(track);
        track_id_counter_eta2++;

        /*
        std::cout << "******************************************************************" << std::endl;
        std::cout << "Resulting track: " << std::endl;
        for (Hit* track_hit : track.getHits()) {
            std::cout << track_hit->getIdx() << ": layer = " << track_hit->getLayer() << ", strip = " << track_hit->getStrip() << ", time = " << track_hit->getTimeEta2() << std::endl;
        }
        std::cout << "\n";
        */
    }
}

// Utility function to update efficiency flags
void Event::updateEfficiencyFlags(const int dt_max, const int dt_min) {

    // Iterate over hits to set efficiency flags
    for (Hit& hit : _hits) {

        // Side η1 efficiency flag updates
        if (hit.hasEta1Time() && hit.getChannel() != _trigger_channel && hit.getRise() == 1) {

            bool within_time_window = true;
            if (_use_external_trigger) {
                if (_trigger_time == -1) {
                    continue;
                }
                within_time_window = (hit.getTimeEta1() - _trigger_time > dt_min && hit.getTimeEta1() - _trigger_time < dt_max);
            } else {
                // If not using external trigger, just check if hit time is within the time window for efficiency consideration
                within_time_window = (hit.getTimeEta1() > dt_min && hit.getTimeEta1() < dt_max);
            }

            // Check if the hit is within the time window for efficiency consideration
            if (within_time_window) {
                _efficiency_flags.eta1_layer[hit.getLayer()] = true;
            } else {
                continue; // Skip hits that are outside the time window for efficiency consideration
            }

            // Additionally check if the hit is part of a valid track on the η1 side to set the track efficiency flag for the layer
            int track_id = hit.getTrackIDEta1();
            auto it = std::find_if(_tracks_eta1.begin(), _tracks_eta1.end(), [track_id](const Track& track) {
                return track.getTrackID() == track_id;
            });

            if (it != _tracks_eta1.end() && it->isValidTrack()) {
                _efficiency_flags.eta1_layer_track[hit.getLayer()] = true;
            }
        }

        // Side η2 efficiency flag updates (same logic as for η1)
        if (hit.hasEta2Time() && hit.getChannel() != _trigger_channel && hit.getRise() == 1) {
            bool within_time_window = true;
            if (_use_external_trigger) {
                if (_trigger_time == -1) {
                    continue;
                }
                within_time_window = (hit.getTimeEta2() - _trigger_time > dt_min && hit.getTimeEta2() - _trigger_time < dt_max);
            } else {
                within_time_window = (hit.getTimeEta2() > dt_min && hit.getTimeEta2() < dt_max);
            }
            if (within_time_window) {
                _efficiency_flags.eta2_layer[hit.getLayer()] = true;
            }

            int track_id = hit.getTrackIDEta2();
            auto it = std::find_if(_tracks_eta2.begin(), _tracks_eta2.end(), [track_id](const Track& track) {
                return track.getTrackID() == track_id;
            });

            if (it != _tracks_eta2.end() && it->isValidTrack()) {
                _efficiency_flags.eta2_layer_track[hit.getLayer()] = true;
            }
        }
    }
}

// Utility function to update efficiency counters based on the efficiency flags set for the event
void Event::updateEfficiencyCounters() {

    // Skip events with no valid trigger time information
    if (_use_external_trigger && _trigger_time == -1) return;

    // External trigger efficiency counter update
    _efficiency_counters.triggered_events_external++;
    _efficiency_counters_tracks.track_triggered_events_external++;
    for (int layer = 0; layer < 3; layer++) {
        if (_efficiency_flags.eta1_layer[layer]) _efficiency_counters.eta1_efficiency_counter[layer]++;
        if (_efficiency_flags.eta2_layer[layer]) _efficiency_counters.eta2_efficiency_counter[layer]++;
        if (_efficiency_flags.eta1_layer[layer] || _efficiency_flags.eta2_layer[layer]) _efficiency_counters.eta_or_efficiency_counter[layer]++;
        if (_efficiency_flags.eta1_layer[layer] && _efficiency_flags.eta2_layer[layer]) _efficiency_counters.eta_and_efficiency_counter[layer]++;

        if (_efficiency_flags.eta1_layer_track[layer]) _efficiency_counters_tracks.track_eta1_efficiency_counter[layer]++;
        if (_efficiency_flags.eta2_layer_track[layer]) _efficiency_counters_tracks.track_eta2_efficiency_counter[layer]++;
        if (_efficiency_flags.eta1_layer_track[layer] || _efficiency_flags.eta2_layer_track[layer]) _efficiency_counters_tracks.track_eta_or_efficiency_counter[layer]++;
        if (_efficiency_flags.eta1_layer_track[layer] && _efficiency_flags.eta2_layer_track[layer]) _efficiency_counters_tracks.track_eta_and_efficiency_counter[layer]++;
    }

    // External trigger + RPC as a trigger efficiency counting
    for (int layer = 0; layer < 3; layer++) {
        const int reference_layer1 = (layer + 1) % 3;
        const int reference_layer2 = (layer + 2) % 3;

        // Check for RPC trigger condition
        const bool rpc_triggered =
            (_efficiency_flags.eta1_layer[reference_layer1] && _efficiency_flags.eta1_layer[reference_layer2]) ||
            (_efficiency_flags.eta2_layer[reference_layer1] && _efficiency_flags.eta2_layer[reference_layer2]);
        const bool rpc_triggered_track =
            (_efficiency_flags.eta1_layer_track[reference_layer1] && _efficiency_flags.eta1_layer_track[reference_layer2]) ||
            (_efficiency_flags.eta2_layer_track[reference_layer1] && _efficiency_flags.eta2_layer_track[reference_layer2]);

        // Update RPC trigger efficiency counters
        if (rpc_triggered) {
            _efficiency_counters.triggered_events_rpc[layer]++;
            if (_efficiency_flags.eta1_layer[layer]) _efficiency_counters.eta1_efficiency_counter_rpc[layer]++;
            if (_efficiency_flags.eta2_layer[layer]) _efficiency_counters.eta2_efficiency_counter_rpc[layer]++;
            if (_efficiency_flags.eta1_layer[layer] || _efficiency_flags.eta2_layer[layer]) _efficiency_counters.eta_or_efficiency_counter_rpc[layer]++;
            if (_efficiency_flags.eta1_layer[layer] && _efficiency_flags.eta2_layer[layer]) _efficiency_counters.eta_and_efficiency_counter_rpc[layer]++;
        }
        if (rpc_triggered_track) {
            _efficiency_counters_tracks.track_triggered_events_rpc[layer]++;
            if (_efficiency_flags.eta1_layer_track[layer]) _efficiency_counters_tracks.track_eta1_efficiency_counter_rpc[layer]++;
            if (_efficiency_flags.eta2_layer_track[layer]) _efficiency_counters_tracks.track_eta2_efficiency_counter_rpc[layer]++;
            if (_efficiency_flags.eta1_layer_track[layer] || _efficiency_flags.eta2_layer_track[layer]) _efficiency_counters_tracks.track_eta_or_efficiency_counter_rpc[layer]++;
            if (_efficiency_flags.eta1_layer_track[layer] && _efficiency_flags.eta2_layer_track[layer]) _efficiency_counters_tracks.track_eta_and_efficiency_counter_rpc[layer]++;
        }
    }
}