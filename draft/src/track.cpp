#include "track.hpp"


// ============================================================
// Track class implementation
// ============================================================

/// Constructor for Track
Track::Track(int id, Hit* first_hit)
    : track_id(id) {
    if (first_hit != nullptr) {
        hits.push_back(first_hit);
        // Initialize layer occupancy based on first hit
        int layer = first_hit->getLayer();
        // TODO: Implement proper layer determination for eta1 and eta2
    }
}

/// Add a hit to the track if it aligns with existing hits
bool Track::addHit(Hit* hit) {
    if (hit == nullptr) return false;
    
    // TODO: Implement hit alignment logic
    // For now, accept all hits
    hits.push_back(hit);
    return true;
}

/// Get whether track has a specific eta1 layer
bool Track::hasEta1Layer(int layer) const {
    if (layer < 0 || layer >= LAYER_COUNT) return false;
    return eta1_layers[layer];
}

/// Get whether track has a specific eta2 layer
bool Track::hasEta2Layer(int layer) const {
    if (layer < 0 || layer >= LAYER_COUNT) return false;
    return eta2_layers[layer];
}

/// Check if track has all eta1 layers
bool Track::hasAllEta1Layers() const {
    for (int i = 0; i < LAYER_COUNT; ++i) {
        if (!eta1_layers[i]) return false;
    }
    return true;
}

/// Check if track has all eta2 layers
bool Track::hasAllEta2Layers() const {
    for (int i = 0; i < LAYER_COUNT; ++i) {
        if (!eta2_layers[i]) return false;
    }
    return true;
}

/// Check if track has all layers for both coordinates
bool Track::hasAllLayers() const {
    return hasAllEta1Layers() && hasAllEta2Layers();
}

/// Get the layer count for a specific coordinate
int Track::getLayerCount(int coordinate) const {
    int count = 0;
    if (coordinate == 0) {  // eta1
        for (int i = 0; i < LAYER_COUNT; ++i) {
            if (eta1_layers[i]) count++;
        }
    } else if (coordinate == 1) {  // eta2
        for (int i = 0; i < LAYER_COUNT; ++i) {
            if (eta2_layers[i]) count++;
        }
    }
    return count;
}
