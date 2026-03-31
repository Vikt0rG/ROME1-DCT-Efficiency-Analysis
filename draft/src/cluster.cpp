#include "cluster.hpp"


// ============================================================
// Cluster class implementation
// ============================================================
/// Constructor that initializes the cluster with the first hit
Cluster::Cluster(Hit* first_hit) : cluster_id(0), center_hit_index(0), tot1(-1), tot2(-1) {
    hits.push_back(first_hit);
}

/// WIP: Core clustering logic: Add a hit to the cluster if it is within the time window and strip distance
bool Cluster::addHit(Hit* hit) {
    bool found_within_strip_distance = false;
    bool found_within_time_window = false;

    // Hits loop to loop for potential cluster partners
    for (size_t idx = 0; idx < hits.size(); idx++) {
        // Check if hit is on the same or neighboring strip
        int strip_diff = abs(hit->getStrip() - hits[idx]->getStrip());
        if (strip_diff <= MAX_STRIP_DISTANCE) {
            found_within_strip_distance = true;
        }
        
        // Check if hit is within time window
        int time_diff = abs(hit->getTimeEta1() - hits[idx]->getTimeEta1());
        if (time_diff < TIME_WINDOW_TICKS) {
            found_within_time_window = true;
            
            // Update cluster center if new hit has smaller rising time
            if (hit->getTimeEta1() < hits[center_hit_index]->getTimeEta1()) {
                center_hit_index = idx;
            }
        }
        
        // Break early if both conditions are satisfied
        if (found_within_strip_distance && found_within_time_window) {
            break;
        }
    }
    
    // Only add hit if both checks pass
    if (found_within_strip_distance && found_within_time_window) {
        hits.push_back(hit);
        return true;
    }
    
    return false;
}