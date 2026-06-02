#include "hit.hpp"


// ==========================================================================================
// Hit class implementation
// ==========================================================================================
/// Constructor that takes raw DCT word and BC0 for processing
Hit::Hit(int clk, int channel, int raw_bcid, int raw_time_eta1, int raw_time_eta2, int rise) 
    // Raw values from parameters
    : _clk(clk)
    , _channel(channel)
    , _raw_bcid(raw_bcid)
    , _raw_time_eta1(raw_time_eta1)
    , _raw_time_eta2(raw_time_eta2)
    , _rise(rise)
{
    geometryMapping();
}

/*
int mod256(int x) { return (x % 256 + 256) % 256; }
    
/// Utility function for applying BC wrap-around correction to BCID
void Hit::applyBCWrapAround(int BC0) {
    _bcid = mod256(_raw_bcid - BC0);
    if (_bcid < -128) _bcid += 256;
    if (_bcid > 128) _bcid -= 256;
}
*/

/// Utility function for mapping channel number to detector layer and strip
void Hit::geometryMapping() {
    int layer = (_channel % 24) / 8;         // Detector layer (0, 1 or 2)
    int column = _channel / 24;              // Detector column (0, 1, 2, ... 5)
    int strip = 8 * column + _channel % 8;   // Strip number (0-47)

    this->_layer = layer;
    this->_column = column;
    this->_strip = strip;
}

/// NEW: WIP: Apply BCID correction after setting BC0 to calculate physical time
void Hit::applyBCIDCorrection(int BC0) {
    // Apply BCID wrap-around
    _bcid = (_raw_bcid - BC0) % 256;
    if (_bcid < -128) _bcid += 256;
    if (_bcid > 128) _bcid -= 256;

    // Convert raw to physical time
    _time_eta1 = TimeUtils::convertRawTimeToPhysical(_raw_time_eta1, _bcid);
    _time_eta2 = TimeUtils::convertRawTimeToPhysical(_raw_time_eta2, _bcid);
} 