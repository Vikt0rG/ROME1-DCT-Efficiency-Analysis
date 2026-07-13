#pragma once

#include <vector>
#include <array>
#include <tuple>

#include "core/types.hpp"
#include "core/hit.hpp"
#include "core/constants.hpp"

// ==========================================================================================
// Track Class: Hits aligned across detector layers
// ==========================================================================================
class Track {
public:
    enum EtaSide { ETA1 = 0, ETA2 = 1 };

    // Constructor
    Track(Hit* first_hit, EtaSide side);

    // Track reconstruction logic
    bool addHit(Hit* hit);
    bool isWithinStripWindow(Hit* hit) const;
    double projectStripPosition(int layer) const;

    // Accessors
    int getTrackID() const { return _track_id; }
    EtaSide getSide() const { return _eta_side; }
    const std::vector<Hit*>& getHits() const { return _track_hits; }

    // Setters
    void setTrackID(int id) { _track_id = id; }

    // Layer occupancy queries (for efficiency calculations)
    bool hasEta1Layer(int layer) const;
    bool hasEta2Layer(int layer) const;
    bool hasAllEta1Layers() const;
    bool hasAllEta2Layers() const;
    bool hasAllLayers() const;  // Both eta1 and eta2 in all layers

    // Track quality methods
    bool isValidTrack() const { return _track_hits.size() >= MIN_TRACK_LENGTH; }
    int getNHits() const;
    int getWidth() const;
    int getLayerCount() const;

    // Track timing information
    int getTimeSeparation() const;  // Time difference between earliest and latest hit in track for the track's eta side
    std::array<std::tuple<bool, bool, int>, LAYER_COUNT> getDts() const;  // Time differences between different layers for the track's eta side
    std::vector<int> getTimeResolution() const;  // Time resolutions for each layer pair

private:
    int _track_id;                  // Unique identifier for the track
    std::vector<Hit*> _track_hits;   // All hits in track (NOT copied, referenced)
    EtaSide _eta_side;               // Track which eta side this track belongs to (ETA1 or ETA2)

    // Layer occupancy flags (for efficiency calculations)
    bool _eta1_layers[3] = {false, false, false};    // Layer 0, 1, 2 for eta1 side
    bool _eta2_layers[3] = {false, false, false};    // Layer 0, 1, 2 for eta2 side
};