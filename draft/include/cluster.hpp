#pragma once

#include "types.hpp"
#include "hit.hpp"


// ============================================================
// Cluster Class: Group of consecutive hits in space/time
// ============================================================
class Cluster {
private:
    int cluster_id;               // Unique identifier for the cluster
    std::vector<Hit*> hits;       // Hits in this cluster (NOT copied, referenced)
    int center_hit_index;         // Index of first-in-time hit (cluster center)
    int tot1, tot2;               // Time-over-threshold for each coordinate

    static constexpr int TIME_WINDOW_TICKS = 18;  // 15 ns clustering time window
    static constexpr int MAX_STRIP_DISTANCE = 1;  // Only consecutive strips

public:
    // Constructor
    Cluster(Hit* first_hit);

    // Core clustering logic
    bool addHit(Hit* hit);

    // Accessors
    Hit* getCenterHit() const { return hits[center_hit_index]; }
    int getSize() const { return hits.size(); }
    int getStrip() const { return getCenterHit()->getStrip(); }
    int getLayer() const { return getCenterHit()->getLayer(); }
    int getCenterTime() const { return getCenterHit()->getTimeEta1(); }
    int getTot1() const { return tot1; }
    int getTot2() const { return tot2; }

    const std::vector<Hit*>& getHits() const { return hits; }

    // ToT calculation (pair rising/falling edges)
    void calculateToT();
};