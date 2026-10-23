#pragma once

#include <vector>
#include <tuple>
#include <array>
#include <string>

struct ColumnsShift {
    int start;
    int end;
    int shift;
};

// =============================================================================
// RUN ENVIRONMENT SWITCH: Set to 1 for TEST_BEAM or 0 for COSMICS
// =============================================================================
#define TEST_BEAM 0
#define BIS2 1

#if TEST_BEAM
    // -------------------------------------------------------------------------
    // TEST BEAM CONFIGURATION
    // ETA_INVERTED = false means that eta1 corresponds to eta+ side
    // -------------------------------------------------------------------------
    inline const std::vector<std::tuple<int, int, int>> DEAD_STRIPS = {{0, 7, 1}};
    constexpr std::array<ColumnsShift, 1> COLUMN_SHIFTS = {{{16, 31, -8}}};

    constexpr bool USE_EXTERNAL_TRIGGER = false;
    constexpr int N_BINS_DT             = 34;

    constexpr int STRIPS_PER_LAYER      = 24;
    constexpr double DETECTOR_LENGTH_CM = 240.0;

    constexpr int MIN_TRACK_LENGTH      = 2;
    constexpr int MAX_STRIP_WINDOW      = 2;
    constexpr int MAX_TIME_WINDOW       = 5;

    constexpr int DEFAULT_DT_MIN        = 0;
    constexpr int DEFAULT_DT_MAX        = 20;

    constexpr bool IS_COSMICS_RUN       = false;
    constexpr bool ETA_INVERTED         = false;

#else
    // -------------------------------------------------------------------------
    // COSMICS CONFIGURATION
    // ETA_INVERTED = true means that eta1 corresponds to eta- side
    // -------------------------------------------------------------------------
        constexpr bool USE_EXTERNAL_TRIGGER = true;
    constexpr int N_BINS_DT             = 22;

    #if BIS2
        constexpr int STRIPS_PER_LAYER  = 40;
        constexpr std::array<ColumnsShift, 1> COLUMN_SHIFTS = {{8, 47, -8}};
        inline const std::vector<std::tuple<int, int, int>> DEAD_STRIPS = {};
    #else
        constexpr int STRIPS_PER_LAYER  = 48;
        constexpr std::array<ColumnsShift, 1> COLUMN_SHIFTS = {};
        inline const std::vector<std::tuple<int, int, int>> DEAD_STRIPS = {};
    #endif
    constexpr double DETECTOR_LENGTH_CM = 160.0;

    constexpr int MIN_TRACK_LENGTH      = 2;
    constexpr int MAX_STRIP_WINDOW      = 2;
    constexpr int MAX_TIME_WINDOW       = 7;

    constexpr int DEFAULT_DT_MIN        = -200;
    constexpr int DEFAULT_DT_MAX        = -100;

    constexpr bool IS_COSMICS_RUN       = true;
    constexpr bool ETA_INVERTED         = true;
#endif

// =============================================================================
// GLOBAL DETECTOR & READOUT CONSTANTS
// =============================================================================
constexpr double SIGNAL_SPEED        = 17.0;    // Signal propagation speed in cm/ns
constexpr int LAYER_COUNT            = 3;       // Number of detector layers (L0, L1, L2)
constexpr int LAYER_PAIR_COUNT       = 3;       // Number of unique layer pairs (L0-L1, L0-L2, L1-L2)
constexpr double STRIP_WIDTH_CM      = 2.5;     // Width of each strip in cm

// Clusterization constants
constexpr int CLUSTER_TIME_WINDOW    = 18;      // Time window for cluster formation (18 ticks = 15 ns)
constexpr int MAX_STRIP_DISTANCE     = 1;       // Maximum strip separation in a cluster

// DCT protocol constants
constexpr int TRIGGER_CHANNEL        = 143;
constexpr int EMPTY_WORD             = 0x5555555;
constexpr int BC_FRAME_SIZE          = 128;
constexpr double TRIGGER_TIME_WINDOW = 400e-9;  // 400 ns
constexpr double TIME_TICK_NS        = 0.833;   // Time resolution of the DCT in ns

// Physical constraints and bounds
constexpr double TOT_ROI_MIN         = 4.0;
constexpr double TOT_ROI_MAX         = 20.0;
constexpr double TR_FIT_CUTOFF_EFF   = 0.25;    // Minimum efficiency threshold for valid TR fit

// Anchor/Partner/Probe approach for efficiency map
constexpr int STRIP_WINDOW          = 1;        // Maximum strip deviation of the anchor w.r.t. the probe
constexpr int DT_WINDOW             = 1;        // Maximum time deviation of the anchor w.r.t. the probe [Ticks]

// Suffixes and labels for plotting and analysis
inline const std::array<std::string, LAYER_PAIR_COUNT> LAYER_PAIR_SUFFIXES = {"0_1", "0_2", "1_2"};