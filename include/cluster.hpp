#pragma once

#include <iostream>
#include "types.hpp"
#include "hit.hpp"
#include "constants.hpp"


// ============================================================
// Cluster Class: Group of consecutive hits in space/time
// ============================================================
class Cluster {
public:
    enum EtaSide { ETA1 = 0, ETA2 = 1 };

private:
    int _cluster_id = -1;                  // Unique identifier for the cluster (set to -1 as invalid amd later updated during clusterization)
    int _layer = -1;                       // Layer which the cluster belongs to (0, 1, or 2) (set to -1 as invalid default)
    std::vector<Hit*> cluster_hits;        // Hits in this cluster (NOT copied, referenced)
    int _center_hit_index = -1;            // Index of first-in-time hit (cluster center) (initialized to -1 as invalid default and updated during clustering)
    int _tot1 = -1, _tot2 = -1;            // Time-over-threshold for each eta side (initialized to -1 as invalid default and updated during ToT calculation)
    EtaSide _eta_side;                     // Track which eta side this cluster belongs to (ETA1 or ETA2)

public:
    // Constructor
    Cluster(Hit* first_hit, EtaSide side, int layer);

    // Core clustering logic
    bool addHit(Hit* hit);

    // ToT calculation for cluster centers
    void calculateToTCluster();

    // Accessors
    Hit* getCenterHit() const { return cluster_hits[_center_hit_index]; }
    EtaSide getSide() const { return _eta_side; }
    int getSize() const { return cluster_hits.size(); }
    int getStrip() const { return getCenterHit()->getStrip(); }
    int getLayer() const { return _layer; }
    int getCenterTime() const { return getCenterHit()->getTimeEta1(); }
    int getTot1() const { return _tot1; }
    int getTot2() const { return _tot2; }
    void setTot1(int value) { _tot1 = value; }
    void setTot2(int value) { _tot2 = value; }

    const std::vector<Hit*>& getHits() const { return cluster_hits; }
};