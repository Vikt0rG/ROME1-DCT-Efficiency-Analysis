#pragma once

#include "core/types.hpp"

class TFile;
class TTree;

/// @namespace summaryHelpers
/// @brief Namespace for helper functions to calculate summary statistics for each
/// measurement entry
namespace summaryHelpers {

    /// @struct Accumulator
    /// @brief Struct to accumulate sums and counts for calculating averages and errors
    struct Accumulator {
        double sum{0.0};  ///< Sum of values accumulated
        size_t hits{0};   ///< Count of hits accumulated
    };

    void setupBranches(TTree* summary_tree, MeasurementMetadata& metadata, MeasurementData& data);

    void getEfficiency(TFile* input_file,
        EfficiencyResults& efficiency_results, EfficiencyResults& efficiency_results_tracks);
    void getClusterSize(TFile* input_file, ClusterSizeResults& cluster_size_results);
    void getRate(TFile* input_file, RateResults& rate_results);
    void getDeadStrips(DeadStrips& dead_strips);
    void getAverageToT(TFile* input_file, ToTResults& tot_results, bool in_beam_only,
        bool in_valid_track_only, DeadStrips& dead_strips, const RoI& region_of_interest);
    void getAverageMultiplicity(TFile* input_file, MultiplicityResults& mult_results,
        bool in_valid_track_only, DeadStrips& dead_strips);
    void processToF(TFile* input_file, ToFResults& tof_results,
        TimeResolutionResults& time_resolution_results, DeadStrips& dead_strips);
}