#pragma once

#include "utils.hpp"

// ============================================================
// Hit Class: Represents a single detector hit
// ============================================================
class Hit {
private:
    // Raw data from the DCT
    int _clk, _channel, _raw_bcid;
    int _raw_time_eta1, _raw_time_eta2, _rise;

    // Processed data - initialized to invalid values until calculated
    int _bcid = -1;                             // After BC wrap-around correction
    int _layer = -1, _column = -1, _strip = -1; // Geometry info
    int _time_eta1 = -1, _time_eta2 = -1;       // Converted times
    int _tot1 = -1, _tot2 = -1;                 // Time-over-threshold

    // IDs - initialized to invalid indices until assigned
    int _hit_id = -1;                   // Unique identifier for this hit
    int _cluster_id_eta1 = -1;          // ID of cluster this hit belongs to (if any)
    int _cluster_id_eta2 = -1;
    int _track_id_eta1 = -1;            // ID of track this hit belongs to (if any)
    int _track_id_eta2 = -1;

public:
    // Constructor from raw word
    Hit(int clk, int channel, int raw_bcid, int raw_time_eta1, int raw_time_eta2, int rise);

    // Getters for hit IDs
    int getIdx() const { return _hit_id; }
    int getClusterIDEta1() const { return _cluster_id_eta1; }
    int getClusterIDEta2() const { return _cluster_id_eta2; }
    int getTrackIDEta1() const { return _track_id_eta1; }
    int getTrackIDEta2() const { return _track_id_eta2; }

    // Getters for raw hit data
    int getClk() const { return _clk; }
    int getChannel() const { return _channel; }
    int getRawBCID() const { return _raw_bcid; }
    int getRawTimeEta1() const { return _raw_time_eta1; }
    int getRawTimeEta2() const { return _raw_time_eta2; }
    int getRise() const { return _rise; }

    // Getters for processed hit data
    int getLayer() const { return _layer; }
    int getStrip() const { return _strip; }
    int getBCID() const { return _bcid; }
    int getTimeEta1() const { return _time_eta1; }
    int getTimeEta2() const { return _time_eta2; }
    int getToT1() const { return _tot1; }
    int getToT2() const { return _tot2; }

    // Setters
    void setIdx(int idx) { _hit_id = idx; }
    void setTot1(int tot) { _tot1 = tot; }
    void setTot2(int tot) { _tot2 = tot; }
    void setClusterIDEta1(int cluster_id) { this->_cluster_id_eta1 = cluster_id; }
    void setClusterIDEta2(int cluster_id) { this->_cluster_id_eta2 = cluster_id; }
    void setTrackIDEta1(int track_id) { this->_track_id_eta1 = track_id; }
    void setTrackIDEta2(int track_id) { this->_track_id_eta2 = track_id; }

    // Processing methods
    void geometryMapping();
    void applyBCIDCorrection(int BC0);

    // Predicates for clustering and track reconstruction
    bool inClusterEta1() const { return _cluster_id_eta1 != -1; }
    bool inClusterEta2() const { return _cluster_id_eta2 != -1; }
    bool inTrackEta1() const { return _track_id_eta1 != -1; }
    bool inTrackEta2() const { return _track_id_eta2 != -1; }

    // Query methods
    bool hasEta1Time() const { return _time_eta1 != -1; }
    bool hasEta2Time() const { return _time_eta2 != -1; }
    bool hasTimeInfo() const { return _time_eta1 != -1 || _time_eta2 != -1; }
};