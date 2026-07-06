#include "track.hpp"

#include <climits>


// ==========================================================================================
// Track class implementation
// ==========================================================================================
// Constructor for Track
Track::Track(Hit* first_hit, EtaSide side)
    : _eta_side(side) {

    if (first_hit != nullptr) {
        _track_hits.push_back(first_hit);
        // Mark the layer of the first hit as occupied
        int layer = first_hit->getLayer();
        if (layer >= 0 && layer < LAYER_COUNT) {
            if (side == ETA1 && first_hit->hasEta1Time()) {
                _eta1_layers[layer] = true;
            } else if (side == ETA2 && first_hit->hasEta2Time()) {
                _eta2_layers[layer] = true;
            }
        } else {
            std::cerr << "Warning: First hit has invalid layer " << layer << std::endl;
        }
    }
}

// Add a hit to the track if it aligns with existing hits
bool Track::addHit(Hit* hit) {
    if (hit == nullptr) return false;

    // Get the hit properties
    int hit_layer = hit->getLayer();
    int hit_time = (_eta_side == ETA1) ? hit->getTimeEta1() : hit->getTimeEta2();

    // Check if hit has valid time information for this track's eta side
    if (hit_time == -1) return false;

    // Ignore hits from the same layer as any existing track member
    for (const Hit* track_hit : _track_hits) {
        if (track_hit->getLayer() == hit_layer) return false;
    }

    /// Strip window check with angular window
    if (!isWithinStripWindow(hit)) return false;

    // Find the closest layer already in the track
    Hit* closest_hit = nullptr;
    int min_layer_diff = INT_MAX;
    
    for (Hit* track_hit : _track_hits) {
        int layer_diff = std::abs(hit_layer - track_hit->getLayer());
        if (layer_diff < min_layer_diff) {
            min_layer_diff = layer_diff;
            closest_hit = track_hit;
        }
    }

    // Check time alignment against the closest layer only
    int closest_time = (_eta_side == ETA1) ? closest_hit->getTimeEta1() : closest_hit->getTimeEta2();
    if (std::abs(hit_time - closest_time) <= MAX_TIME_WINDOW) {
        _track_hits.push_back(hit);
        // Mark the layer as occupied on the appropriate eta side
        if (_eta_side == ETA1) {
            _eta1_layers[hit_layer] = true;
        } else if (_eta_side == ETA2) {
            _eta2_layers[hit_layer] = true;
        }
        return true;
    }

    return false;
}

// Check if a hit is within the strip window of the track's projected position
bool Track::isWithinStripWindow(Hit* hit) const {
    if (_track_hits.size() < 2) {
        // Not enough hits to define a trajectory - just check against first hit
        return std::abs(hit->getStrip() - _track_hits[0]->getStrip()) <= MAX_STRIP_WINDOW;
    }
    
    // Project expected strip position based on existing track hits
    double expected_strip = projectStripPosition(hit->getLayer());
    return std::abs(hit->getStrip() - expected_strip) <= MAX_STRIP_WINDOW;
}

// Simple linear extrapolation to project expected strip position at a given layer
double Track::projectStripPosition(int layer) const {
    if (_track_hits.size() < 2) return _track_hits[0]->getStrip();
    
    int layer0 = _track_hits[0]->getLayer();
    int strip0 = _track_hits[0]->getStrip();
    int layer1 = _track_hits[1]->getLayer();
    int strip1 = _track_hits[1]->getStrip();
    
    double slope = (strip1 - strip0) / (double)(layer1 - layer0);
    return strip0 + slope * (layer - layer0);
}

// Get whether track has a specific eta1 layer
bool Track::hasEta1Layer(int layer) const {
    if (layer < 0 || layer >= LAYER_COUNT) return false;
    return _eta1_layers[layer];
}

// Get whether track has a specific eta2 layer
bool Track::hasEta2Layer(int layer) const {
    if (layer < 0 || layer >= LAYER_COUNT) return false;
    return _eta2_layers[layer];
}

// Check if track has all eta1 layers
bool Track::hasAllEta1Layers() const {
    for (int i = 0; i < LAYER_COUNT; ++i) {
        if (!_eta1_layers[i]) return false;
    }
    return true;
}

// Check if track has all eta2 layers
bool Track::hasAllEta2Layers() const {
    for (int i = 0; i < LAYER_COUNT; ++i) {
        if (!_eta2_layers[i]) return false;
    }
    return true;
}

// Check if track has all layers for both coordinates
bool Track::hasAllLayers() const {
    return hasAllEta1Layers() && hasAllEta2Layers();
}

// Get number of hits in the track
int Track::getNHits() const {
    int count = 0;
    for (const Hit* hit : _track_hits) {
        if (_eta_side == ETA1 && hit->hasEta1Time()) count++;
        else if (_eta_side == ETA2 && hit->hasEta2Time()) count++;
    }
    return count;
}

// Get track width (strip range)
int Track::getWidth() const {
    int min_strip = INT_MAX;
    int max_strip = INT_MIN;
    
    for (const Hit* hit : _track_hits) {
        // Only consider hits with valid times
        if ((_eta_side == ETA1 && !hit->hasEta1Time()) || (_eta_side == ETA2 && !hit->hasEta2Time())) {
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

// Get layer count
int Track::getLayerCount() const {
    int count = 0;
    if (_eta_side == ETA1) {
        for (int i = 0; i < LAYER_COUNT; ++i) {
            if (_eta1_layers[i]) count++;
        }
    } else if (_eta_side == ETA2) {
        for (int i = 0; i < LAYER_COUNT; ++i) {
            if (_eta2_layers[i]) count++;
        }
    }
    return count;
}

// Get track's time differences between layers for the track's eta side
std::array<std::tuple<bool, bool, int>, LAYER_COUNT> Track::getDts() const {
    std::map<int, int> layer_times;

    for (const Hit* hit : _track_hits) {
        int hit_time = (_eta_side == ETA1) ? hit->getTimeEta1() : hit->getTimeEta2();
        if (hit_time == -1) {
            std::cerr << "Warning: Hit in track has invalid time." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        int hit_layer = hit->getLayer();
        if (layer_times.find(hit_layer) == layer_times.end() || hit_time < layer_times[hit_layer]) {
            layer_times[hit_layer] = hit_time;
        }
    }

    // Initialize result tuples with default values: {has_hits=false, is_adjacent=false, dt=0}
    std::array<std::tuple<bool, bool, int>, LAYER_COUNT> result_tuples;
    result_tuples.fill({false, false, 0});

    // Use a flat tracking counter to safely map combinations to the sequential array slots
    int tuple_idx = 0; 
    
    for (int layer = 0; layer < LAYER_COUNT; ++layer) {
        for (int other_layer = layer + 1; other_layer < LAYER_COUNT; ++other_layer) {

            // Dynamically calculate adjacency (if the gap is exactly 1 layer)
            bool is_adjacent = (other_layer - layer == 1);

            // Dynamically evaluate hits and time difference
            if (layer_times.count(layer) && layer_times.count(other_layer)) {
                int dt = layer_times[layer] - layer_times[other_layer];
                result_tuples[tuple_idx] = {true, is_adjacent, dt};
            } else {
                // If a layer is missing a hit, keep track of adjacency structure but flag has_hits as false
                result_tuples[tuple_idx] = {false, is_adjacent, 0};
            }

            // Move safely to the next array slot
            tuple_idx++;
        }
    }

    return result_tuples;
}

// Get track's time separation (dt between earliest and latest hit) for the track's eta side
int Track::getTimeSeparation() const {
    auto dts = getDts();
    int max_dt = -1; // Fallback default for no valid pairs

    // Unpack the tuple: [bool has_hits, bool is_adjacent, int dt_val]
    for (const auto& [has_hits, _, dt_val] : dts) {
        if (has_hits) { 
            int dt = std::abs(dt_val);
            if (dt > max_dt) {
                max_dt = dt;
            }
        }
    }
    return max_dt;
}

// Get time differences between adjacent layers for the track's eta side
std::vector<int> Track::getTimeResolution() const {
    auto dts = getDts();
    std::vector<int> time_differences;

    // Unpack the tuple elements dynamically
    for (const auto& [has_hits, is_adjacent, dt_val] : dts) {
        // Only collect the delta if layers are physically next to each other AND both recorded a hit
        if (is_adjacent && has_hits) {
            time_differences.push_back(std::abs(dt_val));
        }
    }
    return time_differences;
}