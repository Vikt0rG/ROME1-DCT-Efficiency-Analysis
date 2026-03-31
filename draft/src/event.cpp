#include "event.hpp"


// ============================================================
// Event class implementation
// ============================================================
/// Constructor with move semantics for hits vector
Event::Event(int event_number, std::vector<Hit>&& hits_in) 
    : event_number(event_number), hits(std::move(hits_in)) {
    trigger_time = -1;      // Initialize trigger time to invalid
    trigger_channel = 143;  // Initialize trigger channel
}

/// Destructor for the event class
Event::~Event() {
    // Cleanup code if necessary (currently not needed)
}

/// Utility function to add a hit to extract trigger time
void Event::extractTriggerTime() {
    for (const auto& hit : hits) {
        if (hit.getChannel() == TRIGGER_CHANNEL && hit.getRise() == 1) { // Check for rising edge on trigger channel
            trigger_time = hit.getTimeEta2(); // Trigger signal is on the η2 side
            break; // Assuming only one trigger hit per event, we can break after finding it
        }
    }
}

/// WIP: NEW: Clusterization method to form clusters from hits
void Event::clusterize() {
    // Set cluster ID counter
    int cluster_id_counter_eta1 = 0;

    // Side η1 clusterization
    // Loop over the hits and form clusters based on time and strip adjacency
    for (Hit& hit : hits) {
        // Initialize a new cluster object if the hit is not already in a cluster and has valid time information for η1
        if (hit.getTimeEta1() == -1 || hit.inCluster()) continue;

        Cluster cluster(&hit);
        hit.setClusterID(cluster_id_counter_eta1);
        
        // Loop over the remaining hits to find cluster partners for the current hit
        for (Hit& potential_partner : hits) {
            // Skip the same hit
            if (&potential_partner == &hit) continue;

            // Skip hits without valid time information or already in a cluster
            if (potential_partner.getTimeEta1() == -1 || potential_partner.inCluster()) continue;

            // Check logic for cluster membership (time window and strip adjacency)
            if (cluster.addHit(&potential_partner)) potential_partner.setClusterID(cluster_id_counter_eta1);
        }

        // If no more cluster partners are found for the current hit, increment cluster ID counter for the next cluster
        cluster_id_counter_eta1++;
    }

    // Side η2 clusterization (same logic as for η1)
    int cluster_id_counter_eta2 = 0;
    for (Hit& hit : hits) {
        if (hit.getTimeEta2() == -1 || hit.inCluster()) continue;
        Cluster cluster(&hit);
        hit.setClusterID(cluster_id_counter_eta2);
        
        for (Hit& potential_partner : hits) {
            if (&potential_partner == &hit) continue;
            if (potential_partner.getTimeEta2() == -1 || potential_partner.inCluster()) continue;

            if (cluster.addHit(&potential_partner)) potential_partner.setClusterID(cluster_id_counter_eta2);
        }
        cluster_id_counter_eta2++;
    }
}
