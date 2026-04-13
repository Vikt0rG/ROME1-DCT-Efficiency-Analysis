#include "track.hpp"

#include <climits>


// ============================================================
// Track class implementation
// ============================================================
/// Constructor for Track
Track::Track(int id, Hit* first_hit, EtaSide side)
    : track_id(id), eta_side(side) {
    // Initialize layer occupancy arrays
    for (int i = 0; i < LAYER_COUNT; ++i) {
        eta1_layers[i] = false;
        eta2_layers[i] = false;
    }
    
    if (first_hit != nullptr) {
        track_hits.push_back(first_hit);
        // Mark the layer of the first hit as occupied
        int layer = first_hit->getLayer();
        if (layer >= 0 && layer < LAYER_COUNT) {
            if (side == ETA1 && first_hit->hasEta1Time()) {
                eta1_layers[layer] = true;
            } else if (side == ETA2 && first_hit->hasEta2Time()) {
                eta2_layers[layer] = true;
            }
        }
    }
}

/// Add a hit to the track if it aligns with existing hits
bool Track::addHit(Hit* hit) {
    if (hit == nullptr) return false;
    
    // Get the hit properties
    int hit_layer = hit->getLayer();
    int hit_strip = hit->getStrip();
    int hit_time = (eta_side == ETA1) ? hit->getTimeEta1() : hit->getTimeEta2();
    
    // Check if hit has valid time information for this track's eta side
    if (hit_time == -1) return false;
    
    // Check strip window: must be within 5 strips of the first hit
    int first_strip = track_hits[0]->getStrip();
    if (std::abs(hit_strip - first_strip) > MAX_STRIP_WINDOW) return false;
    
    // Check different layer: hit must come from a different layer than existing track hits
    for (const Hit* track_hit : track_hits) {
        if (track_hit->getLayer() == hit_layer) return false;
    }
    
    // Check time alignment: hit must be within 2 ticks of any existing track member
    for (const Hit* track_hit : track_hits) {
        int track_hit_time = (eta_side == ETA1) ? track_hit->getTimeEta1() : track_hit->getTimeEta2();
        if (std::abs(hit_time - track_hit_time) <= MAX_TIME_WINDOW) {
            // Time alignment satisfied with at least one track member
            track_hits.push_back(hit);
            // Mark the layer as occupied on the appropriate eta side
            if (eta_side == ETA1) {
                eta1_layers[hit_layer] = true;
            } else if (eta_side == ETA2) {
                eta2_layers[hit_layer] = true;
            }
            return true;
        }
    }
    
    return false;
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

/// Get the track size for a specific side
int Track::getSize(int side) const {
    int count = 0;
    for (const Hit* hit : track_hits) {
        if (side == 0 && hit->hasEta1Time()) count++;
        else if (side == 1 && hit->hasEta2Time()) count++;
    }
    return count;
}

/// Get the track width (strip range) for a specific side
int Track::getWidth(int side) const {
    int min_strip = INT_MAX;
    int max_strip = INT_MIN;
    
    for (const Hit* hit : track_hits) {
        // Only consider hits with valid times for the requested side
        if ((side == 0 && !hit->hasEta1Time()) || (side == 1 && !hit->hasEta2Time())) {
            continue;
        }
        
        int strip = hit->getStrip();
        if (strip < min_strip) min_strip = strip;
        if (strip > max_strip) max_strip = strip;
    }
    
    // Return 0 if no valid hits found for this side
    if (min_strip == INT_MAX || max_strip == INT_MIN) return 0;
    
    return max_strip - min_strip + 1;
}

/// Get the layer count for a specific side
int Track::getLayerCount(int side) const {
    int count = 0;
    if (side == 0) {  // eta1
        for (int i = 0; i < LAYER_COUNT; ++i) {
            if (eta1_layers[i]) count++;
        }
    } else if (side == 1) {  // eta2
        for (int i = 0; i < LAYER_COUNT; ++i) {
            if (eta2_layers[i]) count++;
        }
    }
    return count;
}