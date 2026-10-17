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
#define TEST_BEAM 1

#if TEST_BEAM
    // -------------------------------------------------------------------------
    // TEST BEAM CONFIGURATION
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

#else
    // -------------------------------------------------------------------------
    // COSMICS CONFIGURATION
    // -------------------------------------------------------------------------
    inline const std::vector<std::tuple<int, int, int>> DEAD_STRIPS = {};
    constexpr std::array<ColumnsShift, 1> COLUMN_SHIFTS = {};

    constexpr bool USE_EXTERNAL_TRIGGER = true;
    constexpr int N_BINS_DT             = 26;

    constexpr int STRIPS_PER_LAYER      = 48;
    constexpr double DETECTOR_LENGTH_CM = 180.0;

    constexpr int MIN_TRACK_LENGTH      = 2;
    constexpr int MAX_STRIP_WINDOW      = 2;
    constexpr int MAX_TIME_WINDOW       = 7;

    constexpr int DEFAULT_DT_MIN        = -200;
    constexpr int DEFAULT_DT_MAX        = -100;

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
constexpr double TRIGGER_TIME_WINDOW = 400e-9; // 400 ns
constexpr double TIME_TICK_NS        = 0.833;  // Time resolution of the DCT in ns

// Physical constraints and bounds
constexpr double TOT_ROI_MIN         = 4.0;
constexpr double TOT_ROI_MAX         = 20.0;

// Suffixes and labels for plotting and analysis
inline const std::array<std::string, LAYER_PAIR_COUNT> LAYER_PAIR_SUFFIXES = {"0_1", "0_2", "1_2"};