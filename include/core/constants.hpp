#pragma once

// ============================================================
// Global Constants for DCT Analysis
// ============================================================

// Detector geometry constants
constexpr int LAYER_COUNT = 3;                      // Number of detector layers (L0, L1, L2)
constexpr int STRIPS_PER_LAYER = 48;                // Number of strips per layer
constexpr float STRIP_WIDTH_CM = 2.5;               // Width of each strip in cm
constexpr float DETECTOR_LENGTH_CM = 180;           // Detector length (cm) for rate calculations

// Clusterization constants
constexpr int CLUSTER_TIME_WINDOW = 18;             // Time window for cluster formation (18 ticks = 15 ns)
constexpr int MAX_STRIP_DISTANCE = 1;               // Maximum strip separation in a cluster (only consecutive strips)

// Track reconstruction constants
constexpr int MIN_TRACK_LENGTH = 2;                 // Minimum number of hits required for a valid track
constexpr int MAX_STRIP_WINDOW = 2;                 // Maximum strip separation for track alignment
constexpr int MAX_TIME_WINDOW = 5;                  // Maximum time difference (ticks) for track alignment

// Time window constants (time difference between hit and trigger signals) for efficiency calculation
constexpr int DEFAULT_DT_MAX = 20;                  // Maximum of the time window (ticks)
constexpr int DEFAULT_DT_MIN = 0;                   // Minimum of the time window (ticks)

// DCT protocol constants
constexpr int TRIGGER_CHANNEL = 143;                // Trigger signal channel
constexpr int EMPTY_WORD = 0x5555555;               // Pattern indicating empty/skipped word
constexpr int BC_FRAME_SIZE = 128;                  // Number of clock cycles per bunch crossing frame (0-127)
constexpr float TRIGGER_TIME_WINDOW = 400e-9;       // Time window for trigger-based rate calculations (400 ns)
constexpr float TIME_TICK_NS = 0.833;               // Duration of one time tick

// Physical contrains and constants
constexpr float TOT_ROI_MIN = 4.0;                  // ToT lower bound for valid hits (ns)
constexpr float TOT_ROI_MAX = 12.0;                 // ToT upper bound for valid hits (ns)