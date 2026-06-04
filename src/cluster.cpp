#include "cluster.hpp"


// ============================================================
// Cluster class implementation
// ============================================================
/// Constructor that initializes the cluster with the first hit
Cluster::Cluster(Hit* first_hit, EtaSide side, int layer)
    : _eta_side(side) {
    cluster_hits.push_back(first_hit);
}

/// Core clustering logic: Add a hit to the cluster if it is within the time window and strip distance
// A hit is added ONLY if it satisfies BOTH strip AND time conditions with AT LEAST ONE existing member
bool Cluster::addHit(Hit* hit) {
    // Check against each existing cluster member
    for (size_t idx = 0; idx < cluster_hits.size(); idx++) {

        // Check if hit is on the same or neighboring strip
        int strip_diff = abs(hit->getStrip() - cluster_hits[idx]->getStrip());
        if (strip_diff > MAX_STRIP_DISTANCE) continue;

        // Check if hit is within time window (use correct eta side)
        int time_diff;
        if (_eta_side == ETA1) {
            time_diff = abs(hit->getTimeEta1() - cluster_hits[idx]->getTimeEta1());
        } else {
            time_diff = abs(hit->getTimeEta2() - cluster_hits[idx]->getTimeEta2());
        }
        if (time_diff >= CLUSTER_TIME_WINDOW) continue;

        // if (eta_side == ETA1) std::cout << "  Potential cluster partner: Hit " << hit->getIdx() << "; Layer: " << hit->getLayer() << "; Strip: " << hit->getStrip() << "; Time " << hit->getTimeEta1() << std::endl;

        // Both conditions satisfied for this cluster member - add hit
        cluster_hits.push_back(hit);

        // Update cluster center if new hit has smaller rising time (use correct eta side)
        if (_eta_side == ETA1) {
            if (hit->getTimeEta1() < cluster_hits[_center_hit_index]->getTimeEta1()) {
                cluster_hits[_center_hit_index]->setIsClusterCenterEta1(false); // Unset old center flag
                _center_hit_index = idx;
                cluster_hits[_center_hit_index]->setIsClusterCenterEta1(true);  // Set new center flag
            }
        } else {
            if (hit->getTimeEta2() < cluster_hits[_center_hit_index]->getTimeEta2()) {
                cluster_hits[_center_hit_index]->setIsClusterCenterEta2(false); // Unset old center flag
                _center_hit_index = idx;
                cluster_hits[_center_hit_index]->setIsClusterCenterEta2(true);  // Set new center flag
            }
        }

        return true;
    }

    return false;
}