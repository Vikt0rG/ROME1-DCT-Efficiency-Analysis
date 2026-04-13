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
    int cluster_id;                       // Unique identifier for the cluster
    std::vector<Hit*> cluster_hits;       // Hits in this cluster (NOT copied, referenced)
    int center_hit_index;                 // Index of first-in-time hit (cluster center)
    int tot1, tot2;                       // Time-over-threshold for both sides
    EtaSide eta_side;                     // Track which eta side this cluster belongs to (ETA1 or ETA2)

public:
    // Constructor
    Cluster(Hit* first_hit, EtaSide side);

    // Core clustering logic
    bool addHit(Hit* hit);

    // ToT calculation for cluster centers
    void calculateToTCluster();

    // Accessors
    Hit* getCenterHit() const { return cluster_hits[center_hit_index]; }
    int getSize() const { return cluster_hits.size(); }
    int getStrip() const { return getCenterHit()->getStrip(); }
    int getLayer() const { return getCenterHit()->getLayer(); }
    int getCenterTime() const { return getCenterHit()->getTimeEta1(); }
    int getTot1() const { return tot1; }
    int getTot2() const { return tot2; }
    void setTot1(int value) { tot1 = value; }
    void setTot2(int value) { tot2 = value; }

    const std::vector<Hit*>& getHits() const { return cluster_hits; }
};