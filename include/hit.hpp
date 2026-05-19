#pragma once

#include "utils.hpp"

// ============================================================
// Hit Class: Represents a single detector hit
// ============================================================
class Hit {
private:
    // Raw data from the DCT
    int clk, channel, raw_bcid, raw_bcout, bcout;
    int raw_time_eta1, raw_time_eta2, rise;

    // Processed data
    int bcid;                     // After BC wrap-around correction
    int layer, column, strip;     // Geometry info
    int time_eta1, time_eta2;     // Converted times
    int tot1, tot2;               // Time-over-threshold

    // IDs
    int hit_id;                   // Unique identifier for this hit
    int cluster_id_eta1;          // ID of cluster this hit belongs to (if any)
    int cluster_id_eta2;
    int track_id_eta1;            // ID of track this hit belongs to (if any)
    int track_id_eta2;

public:
    // Constructor from raw word
    Hit(int clk, int word, int bcout, int BC0);

    // Accessors
    int getIdx() const { return hit_id; }
    int getClusterIDEta1() const { return cluster_id_eta1; }
    int getClusterIDEta2() const { return cluster_id_eta2; }
    int getTrackIDEta1() const { return track_id_eta1; }
    int getTrackIDEta2() const { return track_id_eta2; }

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

    int getToT1() const { return tot1; }
    int getToT2() const { return tot2; }

    // Setters
    void setIdx(int idx) { hit_id = idx; }
    void setTot1(int tot) { tot1 = tot; }
    void setTot2(int tot) { tot2 = tot; }
    void setClusterIDEta1(int cluster_id) { this->cluster_id_eta1 = cluster_id; }
    void setClusterIDEta2(int cluster_id) { this->cluster_id_eta2 = cluster_id; }
    void setTrackIDEta1(int track_id) { this->track_id_eta1 = track_id; }
    void setTrackIDEta2(int track_id) { this->track_id_eta2 = track_id; }

    // Processing methods
    void decodeDCTWord(int word);
    void applyBCWrapAround(int BC0);
    void geometryMapping();       // Maps channel to a layer and strip

    // Predicates for clustering and track reconstruction
    bool inClusterEta1() const { return cluster_id_eta1 != -1; }
    bool inClusterEta2() const { return cluster_id_eta2 != -1; }
    bool inTrackEta1() const { return track_id_eta1 != -1; }
    bool inTrackEta2() const { return track_id_eta2 != -1; }

    // Query methods
    bool hasEta1Time() const { return time_eta1 != -1; }
    bool hasEta2Time() const { return time_eta2 != -1; }
    bool hasTimeInfo() const { return time_eta1 != -1 || time_eta2 != -1; }
};