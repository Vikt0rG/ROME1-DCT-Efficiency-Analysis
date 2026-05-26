#include "hit.hpp"


// ==========================================================================================
// Hit class implementation
// ==========================================================================================
/// Constructor that takes raw DCT word and BC0 for processing
Hit::Hit(int clk, int word, int bcout, int BC0)
    : clk(clk), raw_bcout(bcout) {
    // Decode the raw DCT word to extract channel, BCID, time information and rise/fall edge
    decodeDCTWord(word);

    // Apply BC0 in the event processing
    applyBCWrapAround(BC0);

    // Convert raw to physical time
    time_eta1 = TimeUtils::convertRawTimeToPhysical(raw_time_eta1, bcid);
    time_eta2 = TimeUtils::convertRawTimeToPhysical(raw_time_eta2, bcid);

    // Map channel to layer and strip
    geometryMapping();

    // Initialize ToT values to -1
    tot1 = -1;
    tot2 = -1;

    // Initialize IDs to -1 (indicating not assigned)
    hit_id = -1;
    cluster_id_eta1 = -1;
    cluster_id_eta2 = -1;
    track_id_eta1 = -1;
    track_id_eta2 = -1;
}

/// Utility function for extracting information from the raw DCT word
void Hit::decodeDCTWord(int word) {
    channel = (word >> 20) & 0xFF;          // Bits 20-27
    rise = word & 0x01;                     // Bit 0

    if (rise == 1) {                      // Rising edge
        raw_bcid = word >> 12 & 0xFF;       // Bits 12-19
        raw_time_eta1 = word >> 7 & 0x1F;   // Bits 7-11
        raw_time_eta2 = word >> 1 & 0x1F;   // Bits 1-5

    } else {                              // Falling edge
        raw_bcid = word >> 11 & 0x1FF;      // Bits 11-19
        raw_time_eta1 = word >> 6 & 0x1F;   // Bits 6-10
        raw_time_eta2 = word >> 1 & 0x1F;   // Bits 1-5
    }
}

/// Utility function for applying BC wrap-around correction to BCID and BCOut
void Hit::applyBCWrapAround(int BC0) {
    bcid = raw_bcid - BC0;
    if (bcid < -128) bcid += 256;
    if (bcid > 128) bcid -= 256;

    bcout = (raw_bcout - BC0) % 256;
    if (bcout < -128) bcout += 256;
    if (bcout > 128) bcout -= 256;
}

/*
int mod256(int x) { return (x % 256 + 256) % 256; }

/// Utility function for applying BC wrap-around correction to BCID and BCOut
void Hit::applyBCWrapAround(int BC0) {
    bcid = mod256(raw_bcid - BC0);

    bcout = mod256(raw_bcout - BC0);
    if (bcout < -128) bcout += 256;
    if (bcout > 128) bcout -= 256;
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