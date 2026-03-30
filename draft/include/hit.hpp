#pragma once

#include "utils/utils.hpp"

// ============================================================
// Hit Class: Represents a single detector hit
// ============================================================
class Hit {
private:
    // Raw data from the DCT
    int clk, channel, raw_bcid, raw_bcout, bcout;
    int raw_time_eta1, raw_time_eta2, rise;

    // Processed data
    int bcid, bcout;              // After BC wrap-around correction
    int layer, strip;             // Geometry info
    int time_eta1, time_eta2;     // Converted times

    // IDs (WIP: Decide if neeeded at all)
    int hit_id;                   // Unique identifier for this hit
    int cluster_id;               // ID of cluster this hit belongs to (if any)
    int track_id;                 // ID of track this hit belongs to (if any)

public:
    // Constructor from raw word
    Hit(int clk, int word, int bcout, int BC0);

    // Accessors
    int getClk() const { return clk; }
    int getChannel() const { return channel; }
    int getLayer() const { return layer; }
    int getStrip() const { return strip; }
    int getTimeEta1() const { return time_eta1; }
    int getTimeEta2() const { return time_eta2; }
    int getBCID() const { return bcid; }
    int getBCOut() const { return bcout; }
    int getRawBCID() const { return raw_bcid; }
    int getRawBCOut() const { return raw_bcout; }
    int getRawTimeEta1() const { return raw_time_eta1; }
    int getRawTimeEta2() const { return raw_time_eta2; }
    int getRise() const { return rise; }

    // Processing methods
    void decodeDCTWord(int word);
    void applyBCWrapAround(int BC0);
    void geometryMapping();       // Maps channel to a layer and strip

    // Predicates for clustering and track reconstruction
    bool inCluster();
    bool inTrack();

    // Query methods
    bool hasEta1Time() const { return time_eta1 != -1; }
    bool hasEta2Time() const { return time_eta2 != -1; }
    bool hasTimeInfo() const { return time_eta1 != -1 || time_eta2 != -1; }
};