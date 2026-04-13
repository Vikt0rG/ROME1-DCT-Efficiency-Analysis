#pragma once

#include "types.hpp"
#include "hit.hpp"


// ============================================================
// Track Class: Hits aligned across detector layers
// ============================================================
class Track {
public:
    enum EtaSide { ETA1 = 0, ETA2 = 1 };

private:
    int track_id;                 // Unique identifier for the track
    std::vector<Hit*> track_hits;       // All hits in track (NOT copied, referenced)
    EtaSide eta_side;             // Track which eta side this track belongs to (ETA1 or ETA2)

    // Layer occupancy flags (for efficiency calculations)
    bool eta1_layers[3];          // Layer 0, 1, 2 for η1 coordinate
    bool eta2_layers[3];          // Layer 0, 1, 2 for η2 coordinate

    // Track quality parameters
    static constexpr int MIN_TRACK_LENGTH = 2;    // TODO: Move to constants. Minimum hits required
    static constexpr int LAYER_COUNT = 3;

public:
    // Constructor
    Track(int id, Hit* first_hit, EtaSide side);

    // Track reconstruction logic
    bool addHit(Hit* hit);        // Add hit if aligned with existing hits

    // Accessors
    int getId() const { return track_id; }
    const std::vector<Hit*>& getHits() const { return track_hits; }

    // Layer occupancy queries (for efficiency calculations)
    bool hasEta1Layer(int layer) const;
    bool hasEta2Layer(int layer) const;
    bool hasAllEta1Layers() const;
    bool hasAllEta2Layers() const;
    bool hasAllLayers() const;    // Both eta1 and eta2 in all layers

    // Track quality methods
    bool isValidTrack() const { return track_hits.size() >= MIN_TRACK_LENGTH; }
    int getSize(int side) const;  // side: 0=eta1, 1=eta2
    int getWidth(int side) const;
    int getLayerCount(int side) const;
};