#include "event.hpp"
#include <iostream>
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
        if (hit.getChannel() == trigger_channel && hit.getRise() == 1) { // Check for rising edge on trigger channel
            trigger_time = hit.getTimeEta2(); // Trigger signal is on the η2 side
            break; // Assuming only one trigger hit per event, we can break after finding it
        }
    }
}

/// NEW: Time Over Threshold calculation before clusterization
void Event::calculateTOT() {
    for (Hit& hit : hits) {
        int tot1 = -1;
        int tot2 = -1;
        if (hit.getChannel() == trigger_channel || hit.getRise() == 0) continue; // Skip hits on trigger channel and falling edges
        
        // Find the closest falling edge partner on the same channel (after the rising edge)
        Hit* best_partner_eta1 = nullptr;
        Hit* best_partner_eta2 = nullptr;
        int min_dt_eta1 = INT_MAX;
        int min_dt_eta2 = INT_MAX;
        
        for (Hit& potential_partner : hits) {
            if (&potential_partner == &hit) continue; // Skip the same hit
            if (potential_partner.getChannel() != hit.getChannel() || potential_partner.getRise() != 0) continue;
            
            // Find closest eta1 partner
            if (hit.hasEta1Time() && potential_partner.hasEta1Time()) {
                int dt1 = potential_partner.getTimeEta1() - hit.getTimeEta1();
                if (dt1 > 0 && dt1 < min_dt_eta1) {
                    min_dt_eta1 = dt1;
                    best_partner_eta1 = &potential_partner;
                }
            }
            
            // Find closest eta2 partner
            if (hit.hasEta2Time() && potential_partner.hasEta2Time()) {
                int dt2 = potential_partner.getTimeEta2() - hit.getTimeEta2();
                if (dt2 > 0 && dt2 < min_dt_eta2) {
                    min_dt_eta2 = dt2;
                    best_partner_eta2 = &potential_partner;
                }
            }
        }
        
        if (best_partner_eta1 && min_dt_eta1 > 0) tot1 = min_dt_eta1;
        if (best_partner_eta2 && min_dt_eta2 > 0) tot2 = min_dt_eta2;
        
        hit.setTot1(tot1);
        hit.setTot2(tot2);
    }
}

/// NEW: Clusterization method to form clusters from hits
void Event::clusterize() {
    // Side η1 clusterization
    // Set cluster ID counter
    int cluster_id_counter_eta1 = 0;

    // Loop over the hits and form clusters based on time and strip adjacency
    for (Hit& hit : hits) {

        // Initialize a new cluster object if the rising non-trigger hit is not already in a cluster and has valid time information for η1
        if (hit.getRise() == 0 || hit.getChannel() == trigger_channel || hit.inClusterEta1() || hit.getTimeEta1() == -1) continue;
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "Processing side η1: Hit: " << hit.getIdx() << "; Layer " << hit.getLayer() << "; Strip: " << hit.getStrip() << "; Time " << hit.getTimeEta1() << std::endl;
        Cluster cluster(&hit, Cluster::ETA1);
        hit.setClusterIDEta1(cluster_id_counter_eta1);

        // Loop over the remaining hits to find cluster partners for the current hit
        for (Hit& potential_partner : hits) {
            // Skip the same hit
            if (&potential_partner == &hit) continue;

            // Skip falling hits
            if (potential_partner.getRise() == 0) continue;

            // Skip hits from another layer
            if (potential_partner.getLayer() != hit.getLayer()) continue;

            // Skip hits without valid time information or already in a cluster
            if (potential_partner.getTimeEta1() == -1 || potential_partner.inClusterEta1()) continue;

            // Check logic for cluster membership (time window and strip adjacency)
            if (cluster.addHit(&potential_partner)) potential_partner.setClusterIDEta1(cluster_id_counter_eta1);
        }

        // If no more cluster partners are found for the current hit, store cluster information and increment cluster ID counter for the next cluster
        clusters_eta1.push_back(cluster);
        cluster_id_counter_eta1++;

        // DEBUG: Print out parameters of the hits belonging to the cluster for the current hit
        std::cout << "******************************************************************" << std::endl;
        std::cout << "Resulting cluster: " << std::endl;
        for (Hit* cluster_hit : cluster.getHits()) {
            std::cout << cluster_hit->getIdx() << ": layer = " << cluster_hit->getLayer() << ", strip = " << cluster_hit->getStrip() << ", time = " << cluster_hit->getTimeEta1() << std::endl;
        }
        std::cout << "\n";
    }

    // Side η2 clusterization (same logic as for η1)
    int cluster_id_counter_eta2 = 0;
    for (Hit& hit : hits) {
        if (hit.getRise() == 0 || hit.getChannel() == trigger_channel || hit.getTimeEta2() == -1 || hit.inClusterEta2()) continue;
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "Processing side η2: Hit: " << hit.getIdx() << "; Layer " << hit.getLayer() << "; Strip: " << hit.getStrip() << "; Time " << hit.getTimeEta2() << std::endl;
        Cluster cluster(&hit, Cluster::ETA2);
        hit.setClusterIDEta2(cluster_id_counter_eta2);

        for (Hit& potential_partner : hits) {
            if (&potential_partner == &hit) continue;

            if (potential_partner.getRise() == 0) continue;

            if (potential_partner.getLayer() != hit.getLayer()) continue;

            if (potential_partner.getTimeEta2() == -1 || potential_partner.inClusterEta2()) continue;

            if (cluster.addHit(&potential_partner)) potential_partner.setClusterIDEta2(cluster_id_counter_eta2);
        }
        clusters_eta2.push_back(cluster);
        cluster_id_counter_eta2++;

        // DEBUG: Print out parameters of the hits belonging to the cluster for the current hit
        std::cout << "******************************************************************" << std::endl;
        std::cout << "Resulting cluster: " << std::endl;
        for (Hit* cluster_hit : cluster.getHits()) {
            std::cout << cluster_hit->getIdx() << ": layer = " << cluster_hit->getLayer() << ", strip = " << cluster_hit->getStrip() << ", time = " << cluster_hit->getTimeEta2() << std::endl;
        }
        std::cout << "\n";
    }

    // NOTE: Clusters are retained for later access (not cleared here)
}

/// WIP: Time-over-threshold calculation for each cluster
void Event::calculateTOTCluster() {
    for (Cluster& cluster : clusters_eta1) {
        cluster.calculateToTCluster();
    }
    for (Cluster& cluster : clusters_eta2) {
        cluster.calculateToTCluster();
    }
}

/// Track reconstruction from clusters (WIP)
void Event::reconstructTracks() {
    // TODO: Implement track reconstruction from clusters
    // For now, this is a placeholder to allow compilation
}