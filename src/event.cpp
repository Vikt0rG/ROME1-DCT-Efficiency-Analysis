#include "event.hpp"


// ============================================================
// Event class implementation
// ============================================================
/// Constructor without hits
Event::Event(int event_number, 
    EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks, bool use_external_trigger)
    :   event_number(event_number),
        efficiency_counters(counters),
        efficiency_counters_tracks(counters_tracks),
        use_external_trigger(use_external_trigger) {
    trigger_time = -1;      // Initialize trigger time to invalid
    trigger_channel = TRIGGER_CHANNEL;  // Initialize trigger channel from constants
}

/// Constructor with move semantics for hits vector
Event::Event(int event_number, std::vector<Hit>&& hits_in, 
    EfficiencyCounters& counters, EfficiencyCountersTracks& counters_tracks, bool use_external_trigger) 
    :   event_number(event_number),
        hits(std::move(hits_in)),
        efficiency_counters(counters),
        efficiency_counters_tracks(counters_tracks),
        use_external_trigger(use_external_trigger) {
    trigger_time = -1;      // Initialize trigger time to invalid
    trigger_channel = TRIGGER_CHANNEL;  // Initialize trigger channel from constants

    // Reset efficiency flags for each new event
    for (int i = 0; i < 3; i++) {
        efficiency_flags.eta1_layer[i] = false;
        efficiency_flags.eta2_layer[i] = false;
        efficiency_flags.eta1_layer_track[i] = false;
        efficiency_flags.eta2_layer_track[i] = false;
    }
}

/// Destructor for the event class
Event::~Event() {
    // Cleanup code if necessary (currently not needed)
}

/// Utility function to add a hit to extract trigger time
void Event::extractTriggerTime() {
    if (!use_external_trigger) {
        trigger_time = -1;
        return;
    }
    for (const auto& hit : hits) {
        if (hit.getChannel() == trigger_channel && hit.getRise() == 1) { // Check for rising edge on trigger channel
            trigger_time = hit.getTimeEta2(); // Trigger signal is on the η2 side
            break; // Assuming only one trigger hit per event, we can break after finding it
        }
    }
}

/// Utility function to estimate whether there are more rising or falling edges in the event for ToT calculation
bool hasMoreRisingEdges(const std::vector<Hit>& hits) {
    int rising_count = 0;
    int falling_count = 0;
    for (const auto& hit : hits) {
        if (hit.getRise() == 1) rising_count++;
        else if (hit.getRise() == 0) falling_count++;
    }
    return rising_count > falling_count;
}

/// Time Over Threshold calculation before clusterization
void Event::calculateTOT() {
    // Vectors to keep track which hits have been paired for ToT calculation to avoid double counting
    // Note: For pairing rely on eta1 side only unless it doesn't have time information, then use eta2 side.
    std::vector<bool> paired(hits.size(), false);

    // Determine if we have more rising edges than falling edges in the event, which can indicate potential data issues (?)
    bool more_rising_edges = hasMoreRisingEdges(hits);
    bool more_falling_edges = !more_rising_edges;

    for (size_t i = 0; i < hits.size(); i++) {
        Hit& hit = hits[i];
        int tot1 = -1;
        int tot2 = -1;

        if (hit.getChannel() == trigger_channel) continue; // Skip hits on trigger channel

        if (paired[i]) continue; // Skip already paired hits

        // Greedy approach: Skip hits that are more common, to match those to the less common ones
        if ((hit.getRise() == 0 && more_falling_edges) || (hit.getRise() == 1 && more_rising_edges)) continue;

        // Find the closest in-time partner edge on the same channel
        int best_j = -1;
        int min_dt_eta1 = INT_MAX;
        int min_dt_eta2 = INT_MAX;
        const int target_rise = (hit.getRise() == 1) ? 0 : 1;

        for (size_t j = 0; j < hits.size(); j++) {
            Hit& potential_partner = hits[j];
            if (&potential_partner == &hit || paired[j]) continue;  // Skip the same hit or already paired hits
            if (potential_partner.getChannel() != hit.getChannel()) continue;   // Skip hits from different channels
            if (potential_partner.getRise() != target_rise) continue;   // Skip hits with the same edge type

            // Find closest partner using eta1 time information if available (fallback to eta2 if not)
            if (hit.hasEta1Time() && potential_partner.hasEta1Time()) {
                int dt1 = potential_partner.getTimeEta1() - hit.getTimeEta1();
                if (dt1 > 0 && dt1 < min_dt_eta1) {
                    min_dt_eta1 = dt1;
                    best_j = static_cast<int>(j);
                }
            } else if (hit.hasEta2Time() && potential_partner.hasEta2Time()) {
                int dt2 = potential_partner.getTimeEta2() - hit.getTimeEta2();
                if (dt2 > 0 && dt2 < min_dt_eta2) {
                    min_dt_eta2 = dt2;
                    best_j = static_cast<int>(j);
                }
            }
        }

        if (best_j >= 0) {
            paired[i] = true;
            paired[best_j] = true;
            if (min_dt_eta1 < INT_MAX) tot1 = min_dt_eta1;
            if (min_dt_eta2 < INT_MAX) tot2 = min_dt_eta2;
        }

        hit.setTot1(tot1);
        hit.setTot2(tot2);
    }
}

/// Clusterization method to form clusters from hits
void Event::clusterize() {
    // Side η1 clusterization
    // Set cluster ID counter
    int cluster_id_counter_eta1 = 0;

    // Loop over the hits and form clusters based on time and strip adjacency
    for (Hit& hit : hits) {

        // Initialize a new cluster object if the rising non-trigger hit is not already in a cluster and has valid time information for η1
        if (hit.getRise() == 0 || hit.getChannel() == trigger_channel || hit.inClusterEta1() || hit.getTimeEta1() == -1) continue;
        // std::cout << "------------------------------------------------------------------" << std::endl;
        // std::cout << "Processing side η1: Hit: " << hit.getIdx() << "; Layer " << hit.getLayer() << "; Strip: " << hit.getStrip() << "; Time " << hit.getTimeEta1() << std::endl;
        Cluster cluster(&hit, Cluster::ETA1, hit.getLayer());
        hit.setClusterIDEta1(cluster_id_counter_eta1);

        // Loop over the remaining hits to find cluster partners for the current hit
        for (Hit& potential_partner : hits) {
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
        clusters_eta1.push_back(cluster);
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
    for (Hit& hit : hits) {
        if (hit.getRise() == 0 || hit.getChannel() == trigger_channel || hit.getTimeEta2() == -1 || hit.inClusterEta2()) continue;
        // std::cout << "------------------------------------------------------------------" << std::endl;
        // std::cout << "Processing side η2: Hit: " << hit.getIdx() << "; Layer " << hit.getLayer() << "; Strip: " << hit.getStrip() << "; Time " << hit.getTimeEta2() << std::endl;
        Cluster cluster(&hit, Cluster::ETA2, hit.getLayer());
        hit.setClusterIDEta2(cluster_id_counter_eta2);

        for (Hit& potential_partner : hits) {
            if (&potential_partner == &hit) continue;

            if (potential_partner.getRise() == 0) continue;

            if (potential_partner.getLayer() != hit.getLayer()) continue;

            if (potential_partner.inClusterEta2()) continue;

            if (potential_partner.getTimeEta2() == -1) continue;

            if (cluster.addHit(&potential_partner)) potential_partner.setClusterIDEta2(cluster_id_counter_eta2);
        }
        clusters_eta2.push_back(cluster);
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

/// Time-over-threshold calculation for cluster centers (after clusterization)
void Event::calculateTOTCluster() {
    // Calculate TOT for eta1 cluster centers
    for (Cluster& cluster : clusters_eta1) {
        Hit* center_hit = cluster.getCenterHit();
        int tot1 = -1;
        
        // Find closest falling edge partner on the same channel (after the rising edge)
        Hit* best_partner = nullptr;
        int min_dt = INT_MAX;
        
        for (Hit& potential_partner : hits) {
            if (potential_partner.getChannel() != center_hit->getChannel()) continue;
            if (potential_partner.getRise() != 0) continue;  // Must be falling edge
            
            int dt = potential_partner.getTimeEta1() - center_hit->getTimeEta1();
            if (dt > 0 && dt < min_dt) {
                min_dt = dt;
                best_partner = &potential_partner;
            }
        }
        
        if (best_partner && min_dt > 0) tot1 = min_dt;
        cluster.setTot1(tot1);
    }
    
    // Calculate TOT for eta2 cluster centers
    for (Cluster& cluster : clusters_eta2) {
        Hit* center_hit = cluster.getCenterHit();
        int tot2 = -1;
        
        // Find closest falling edge partner on the same channel (after the rising edge)
        Hit* best_partner = nullptr;
        int min_dt = INT_MAX;
        
        for (Hit& potential_partner : hits) {
            if (potential_partner.getChannel() != center_hit->getChannel()) continue;
            if (potential_partner.getRise() != 0) continue;  // Must be falling edge
            
            int dt = potential_partner.getTimeEta2() - center_hit->getTimeEta2();
            if (dt > 0 && dt < min_dt) {
                min_dt = dt;
                best_partner = &potential_partner;
            }
        }
        
        if (best_partner && min_dt > 0) tot2 = min_dt;
        cluster.setTot2(tot2);
    }
}

/// Track reconstruction from clusters
void Event::reconstructTracks() {
    // Side η1 track reconstruction logic
    int track_id_counter_eta1 = 0;

    // Loop over η1 cluster centers and form tracks based on time and strip alignment across layers
    for (Cluster& cluster : clusters_eta1) {
        // Initialize a new track object if the cluster center is not already in a track
        // NOTE: Here we no longer need checks for rising edge, non-trigger channel and valid
        // time information as these are already ensured for cluster centers at the cluster level
        if (cluster.getCenterHit()->inTrackEta1()) continue;
        // std::cout << "------------------------------------------------------------------" << std::endl;
        // std::cout << "Processing track reconstruction for side η1: " << std::endl;
        Track track(track_id_counter_eta1, cluster.getCenterHit(), Track::ETA1);
        cluster.getCenterHit()->setTrackIDEta1(track_id_counter_eta1);

        // Loop over the remaining cluster centers to find track partners for the current cluster center
        for (Cluster& potential_partner_cluster : clusters_eta1) {
            // Skip the same cluster
            if (&potential_partner_cluster == &cluster) continue;

            // Skip cluster centers from the same layer
            if (potential_partner_cluster.getCenterHit()->getLayer() == cluster.getCenterHit()->getLayer()) continue;

            // Skip cluster centers that are already in a track
            if (potential_partner_cluster.getCenterHit()->inTrackEta1()) continue;

            // Check logic for track membership (strip window and time alignment)
            if (track.addHit(potential_partner_cluster.getCenterHit())) potential_partner_cluster.getCenterHit()->setTrackIDEta1(track_id_counter_eta1);
        }

        // If no more track partners are found for the current cluster center, store track information and increment track ID counter for the next track
        tracks_eta1.push_back(track);
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
    for (Cluster& cluster : clusters_eta2) {
        if (cluster.getCenterHit()->inTrackEta2()) continue;
        // std::cout << "------------------------------------------------------------------" << std::endl;
        // std::cout << "Processing track reconstruction for side η2: " << std::endl;
        Track track(track_id_counter_eta2, cluster.getCenterHit(), Track::ETA2);
        cluster.getCenterHit()->setTrackIDEta2(track_id_counter_eta2);

        for (Cluster& potential_partner_cluster : clusters_eta2) {
            if (&potential_partner_cluster == &cluster) continue;

            if (potential_partner_cluster.getCenterHit()->getLayer() == cluster.getCenterHit()->getLayer()) continue;

            if (potential_partner_cluster.getCenterHit()->inTrackEta2()) continue;

            if (track.addHit(potential_partner_cluster.getCenterHit())) potential_partner_cluster.getCenterHit()->setTrackIDEta2(track_id_counter_eta2);
        }
        tracks_eta2.push_back(track);
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

/// Utility function to update efficiency flags
void Event::updateEfficiencyFlags(const int dt_max, const int dt_min) {

    // Iterate over hits to set efficiency flags
    for (Hit& hit : hits) {

        // Side η1 efficiency flag updates
        if (hit.hasEta1Time() && hit.getChannel() != trigger_channel && hit.getRise() == 1) {

            bool within_time_window = true;
            if (use_external_trigger) {
                if (trigger_time == -1) {
                    continue;
                }
                within_time_window = (hit.getTimeEta1() - trigger_time > dt_min && hit.getTimeEta1() - trigger_time < dt_max);
            } else {
                // If not using external trigger, just check if hit time is within the time window for efficiency consideration
                within_time_window = (hit.getTimeEta1() > dt_min && hit.getTimeEta1() < dt_max);
            }

            // Check if the hit is within the time window for efficiency consideration
            if (within_time_window) {
                efficiency_flags.eta1_layer[hit.getLayer()] = true;
            } else {
                continue; // Skip hits that are outside the time window for efficiency consideration
            }

            // Additionally check if the hit is part of a valid track on the η1 side to set the track efficiency flag for the layer
            int track_id = hit.getTrackIDEta1();
            auto it = std::find_if(tracks_eta1.begin(), tracks_eta1.end(), [track_id](const Track& track) {
                return track.getId() == track_id;
            });

            if (it != tracks_eta1.end() && it->isValidTrack()) {
                efficiency_flags.eta1_layer_track[hit.getLayer()] = true;
            }
        }

        // Side η2 efficiency flag updates (same logic as for η1)
        if (hit.hasEta2Time() && hit.getChannel() != trigger_channel && hit.getRise() == 1) {
            bool within_time_window = true;
            if (use_external_trigger) {
                if (trigger_time == -1) {
                    continue;
                }
                within_time_window = (hit.getTimeEta2() - trigger_time > dt_min && hit.getTimeEta2() - trigger_time < dt_max);
            } else {
                within_time_window = (hit.getTimeEta2() > dt_min && hit.getTimeEta2() < dt_max);
            }
            if (within_time_window) {
                efficiency_flags.eta2_layer[hit.getLayer()] = true;
            }

            int track_id = hit.getTrackIDEta2();
            auto it = std::find_if(tracks_eta2.begin(), tracks_eta2.end(), [track_id](const Track& track) {
                return track.getId() == track_id;
            });

            if (it != tracks_eta2.end() && it->isValidTrack()) {
                efficiency_flags.eta2_layer_track[hit.getLayer()] = true;
            }
        }
    }
}

/// Utility function to update efficiency counters based on the efficiency flags set for the event
void Event::updateEfficiencyCounters() {

    // Skip events with no valid trigger time information
    if (use_external_trigger && trigger_time == -1) return;

    // External trigger efficiency counter update
    efficiency_counters.triggered_events_external++;
    efficiency_counters_tracks.track_triggered_events_external++;
    for (int layer = 0; layer < 3; layer++) {
        if (efficiency_flags.eta1_layer[layer]) efficiency_counters.eta1_efficiency_counter[layer]++;
        if (efficiency_flags.eta2_layer[layer]) efficiency_counters.eta2_efficiency_counter[layer]++;
        if (efficiency_flags.eta1_layer[layer] || efficiency_flags.eta2_layer[layer]) efficiency_counters.eta_or_efficiency_counter[layer]++;
        if (efficiency_flags.eta1_layer[layer] && efficiency_flags.eta2_layer[layer]) efficiency_counters.eta_and_efficiency_counter[layer]++;

        if (efficiency_flags.eta1_layer_track[layer]) efficiency_counters_tracks.track_eta1_efficiency_counter[layer]++;
        if (efficiency_flags.eta2_layer_track[layer]) efficiency_counters_tracks.track_eta2_efficiency_counter[layer]++;
        if (efficiency_flags.eta1_layer_track[layer] || efficiency_flags.eta2_layer_track[layer]) efficiency_counters_tracks.track_eta_or_efficiency_counter[layer]++;
        if (efficiency_flags.eta1_layer_track[layer] && efficiency_flags.eta2_layer_track[layer]) efficiency_counters_tracks.track_eta_and_efficiency_counter[layer]++;
    }

    // External trigger + RPC as a trigger efficiency counting
    for (int layer = 0; layer < 3; layer++) {
        const int reference_layer1 = (layer + 1) % 3;
        const int reference_layer2 = (layer + 2) % 3;
        const bool rpc_triggered =
            (efficiency_flags.eta1_layer[reference_layer1] && efficiency_flags.eta1_layer[reference_layer2]) ||
            (efficiency_flags.eta2_layer[reference_layer1] && efficiency_flags.eta2_layer[reference_layer2]);
        const bool rpc_triggered_track =
            (efficiency_flags.eta1_layer_track[reference_layer1] && efficiency_flags.eta1_layer_track[reference_layer2]) ||
            (efficiency_flags.eta2_layer_track[reference_layer1] && efficiency_flags.eta2_layer_track[reference_layer2]);

        if (rpc_triggered) {
            efficiency_counters.triggered_events_rpc[layer]++;
            if (efficiency_flags.eta1_layer[layer]) efficiency_counters.eta1_efficiency_counter_rpc[layer]++;
            if (efficiency_flags.eta2_layer[layer]) efficiency_counters.eta2_efficiency_counter_rpc[layer]++;
            if (efficiency_flags.eta1_layer[layer] || efficiency_flags.eta2_layer[layer]) efficiency_counters.eta_or_efficiency_counter_rpc[layer]++;
            if (efficiency_flags.eta1_layer[layer] && efficiency_flags.eta2_layer[layer]) efficiency_counters.eta_and_efficiency_counter_rpc[layer]++;
        }

        if (rpc_triggered_track) {
            efficiency_counters_tracks.track_triggered_events_rpc[layer]++;
            if (efficiency_flags.eta1_layer_track[layer]) efficiency_counters_tracks.track_eta1_efficiency_counter_rpc[layer]++;
            if (efficiency_flags.eta2_layer_track[layer]) efficiency_counters_tracks.track_eta2_efficiency_counter_rpc[layer]++;
            if (efficiency_flags.eta1_layer_track[layer] || efficiency_flags.eta2_layer_track[layer]) efficiency_counters_tracks.track_eta_or_efficiency_counter_rpc[layer]++;
            if (efficiency_flags.eta1_layer_track[layer] && efficiency_flags.eta2_layer_track[layer]) efficiency_counters_tracks.track_eta_and_efficiency_counter_rpc[layer]++;
        }
    }
}