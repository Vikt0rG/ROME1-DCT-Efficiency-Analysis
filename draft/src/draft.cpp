#include "draft.hpp"

// Hit class implementation
Hit::Hit(int clk, int word, int raw_bcout) : clk(clk), raw_bcout(raw_bcout) {
    decodeDCTWord(word);
    applyBCWrapAround(0);  // BC0 will be applied later in the event processing
    geometryMapping();
}

void Hit::decodeDCTWord(int word) {
    channel = (word >> 20) & 0xFF;          // Bits 20-27
    rise = word & 0x01;                     // Bit 0

    if (rise == 1) {                      // Rising edge
        raw_bcid = word >> 12 & 0xFF;       // Bits 12-19
        raw_time_eta1 = word >> 7 & 0x1F;   // Bits 7-11
        raw_time_eta2 = word >> 1 & 0x3F;   // Bits 1-5

    } else {                              // Falling edge
        raw_bcid = word >> 11 & 0x1FF;      // Bits 11-19
        raw_time_eta1 = word >> 6 & 0x1F;   // Bits 6-10
        raw_time_eta2 = word >> 1 & 0x1F;   // Bits 1-5
    }
}

void Hit::applyBCWrapAround(int BC0) {
    bcid = raw_bcid - BC0;
    if (bcid < -128) bcid += 256;
    if (bcid > 128) bcid -= 256;

    bcout = raw_bcout - BC0;
    if (bcout < -128) bcout += 256;
    if (bcout > 128) bcout -= 256;
}

void Hit::geometryMapping() {
    int layer = (channel % 24) / 8;         // Detector layer (0, 1 or 2)
    int column = channel / 24;              // Detector column (0, 1, 2, ... 5)
    int strip = 8 * column + channel % 8;   // Strip number (0-47)

    this->layer = layer;
    this->strip = strip;
}

Event::Event(int event_number) : event_number(event_number) {
    trigger_time = -1;      // Initialize trigger time to invalid
    trigger_channel = 143;  // Initialize trigger channel
}
