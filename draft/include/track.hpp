#pragma once

#include "types.hpp"
#include "hit.hpp"


// ============================================================
// Track Class: Hits aligned across detector layers
// ============================================================
class Track {
private:
    std::vector<Hit*> hits;       // All hits in track (NOT copied, referenced)
    int track_id;

    // Layer occupancy flags (for efficiency calculations)
    bool eta1_layers[3];          // Layer 0, 1, 2 for η1 coordinate
    bool eta2_layers[3];          // Layer 0, 1, 2 for η2 coordinate

    // Track quality parameters
    static constexpr int MIN_TRACK_LENGTH = 2;    // Minimum hits required
    static constexpr int LAYER_COUNT = 3;

public:
    // Constructor
    Track(int id, Hit* first_hit);

    // Track reconstruction logic
    bool addHit(Hit* hit);        // Add hit if aligned with existing hits

    // Accessors
    int getId() const { return track_id; }
    int getLength() const { return hits.size(); }
    const std::vector<Hit*>& getHits() const { return hits; }

    // Layer occupancy queries (for efficiency calculations)
    bool hasEta1Layer(int layer) const;
    bool hasEta2Layer(int layer) const;
    bool hasAllEta1Layers() const;
    bool hasAllEta2Layers() const;
    bool hasAllLayers() const;    // Both eta1 and eta2 in all layers

    // Track quality methods
    bool isValidTrack() const { return hits.size() >= MIN_TRACK_LENGTH; }
    int getLayerCount(int coordinate) const;  // coordinate: 0=eta1, 1=eta2
};