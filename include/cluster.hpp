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

    // Constructor
    Cluster(Hit* first_hit, EtaSide side, int layer);

    // Core clustering logic
    bool addHit(Hit* hit);

    // Accessors
    int getClusterID() const { return _cluster_id; }
    EtaSide getSide() const { return _eta_side; }
    Hit* getCenterHit() const { return _center_hit; }
    const std::vector<Hit*>& getHits() const { return _cluster_hits; }
    int getCenterStrip() const { return _center_hit ? _center_hit->getStrip() : -1; }
    int getSize() const { return static_cast<int>(_cluster_hits.size()); }
    int getLayer() const { return _layer; }
    int getCenterTime() const { 
        if (!_center_hit) return -1;
        return (_eta_side == ETA1) ? _center_hit->getTimeEta1() : _center_hit->getTimeEta2();
    }
    int getTot1() const { return _tot1; }
    int getTot2() const { return _tot2; }

    // Setters
    void setClusterID(int id) { _cluster_id = id; }
    void setTot1(int value) { _tot1 = value; }
    void setTot2(int value) { _tot2 = value; }

private:
    int _cluster_id = -1;                  // Unique identifier for the cluster (set to -1 as invalid and later updated during clusterization)
    int _layer = -1;                       // Layer which the cluster belongs to (0, 1, or 2) (set to -1 as invalid default)
    Hit* _center_hit = nullptr;            // Direct pointer to the first-in-time hit in the cluster (cluster center)
    std::vector<Hit*> _cluster_hits;       // Hits in this cluster (NOT copied, referenced)
    int _tot1 = -1, _tot2 = -1;            // Time-over-threshold for each eta side (initialized to -1 as invalid default and updated during ToT calculation)
    EtaSide _eta_side;                     // Track which eta side this cluster belongs to (ETA1 or ETA2)
};