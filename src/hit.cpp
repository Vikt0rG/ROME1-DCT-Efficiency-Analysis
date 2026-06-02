#include "hit.hpp"


// ==========================================================================================
// Hit class implementation
// ==========================================================================================
/// Constructor that takes raw DCT word and BC0 for processing
Hit::Hit(int clk, int channel, int raw_bcid, int raw_time_eta1, int raw_time_eta2, int rise) : clk(clk) {

    // Map channel to layer, column and strip
    geometryMapping();

    // Initialize variables to calculate as -1 (indicating invalid/unavailable) until processed
    bcid = -1;
    time_eta1 = -1;
    time_eta2 = -1;
    tot1 = -1;
    tot2 = -1;

    // Initialize IDs to -1 (indicating not assigned)
    hit_id = -1;
    cluster_id_eta1 = -1;
    cluster_id_eta2 = -1;
    track_id_eta1 = -1;
    track_id_eta2 = -1;
}

/*
int mod256(int x) { return (x % 256 + 256) % 256; }
    
/// Utility function for applying BC wrap-around correction to BCID
void Hit::applyBCWrapAround(int BC0) {
    bcid = mod256(raw_bcid - BC0);
    if (bcid < -128) bcid += 256;
    if (bcid > 128) bcid -= 256;
}
*/

/// Utility function for mapping channel number to detector layer and strip
void Hit::geometryMapping() {
    int layer = (channel % 24) / 8;         // Detector layer (0, 1 or 2)
    int column = channel / 24;              // Detector column (0, 1, 2, ... 5)
    int strip = 8 * column + channel % 8;   // Strip number (0-47)

    this->layer = layer;
    this->column = column;
    this->strip = strip;
}

/// NEW: WIP: Apply BCID correction after setting BC0 to calculate physical time
void Hit::applyBCIDCorrection(int BC0) {
    // Apply BCID wrap-around
    bcid = (raw_bcid - BC0) % 256;
    if (bcid < -128) bcid += 256;
    if (bcid > 128) bcid -= 256;

    // Convert raw to physical time
    time_eta1 = TimeUtils::convertRawTimeToPhysical(raw_time_eta1, bcid);
    time_eta2 = TimeUtils::convertRawTimeToPhysical(raw_time_eta2, bcid);
} 