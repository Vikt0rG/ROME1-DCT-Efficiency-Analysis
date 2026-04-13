#pragma once

// ============================================================
// Global Constants for DCT Analysis
// ============================================================

// Detector geometry constants
constexpr int LAYER_COUNT = 3;              // Number of detector layers (L0, L1, L2)
constexpr int STRIPS_PER_LAYER = 48;        // Number of strips per layer

// Clusterization constants
constexpr int CLUSTER_TIME_WINDOW = 18;     // Time window for cluster formation (18 ticks = 15 ns)
constexpr int MAX_STRIP_DISTANCE = 1;       // Maximum strip separation in a cluster (only consecutive strips)

// Track reconstruction constants
constexpr int MIN_TRACK_LENGTH = 2;         // Minimum number of hits required for a valid track
constexpr int MAX_STRIP_WINDOW = 5;         // Maximum strip separation for track alignment
constexpr int MAX_TIME_WINDOW = 2;          // Maximum time difference (ticks) for track alignment

// Time window constants (time difference between hit and trigger signals) for efficiency calculation
constexpr int DEFAULT_DT_MAX = -100;        // Maximum of the time window (ticks)
constexpr int DEFAULT_DT_MIN = -180;        // Minimum of the time window (ticks)

// DCT protocol constants
constexpr int TRIGGER_CHANNEL = 143;        // Trigger signal channel
constexpr int EMPTY_WORD = 0x5555555;       // Pattern indicating empty/skipped word
constexpr int BC_FRAME_SIZE = 128;          // Number of clock cycles per bunch crossing frame (0-127)
