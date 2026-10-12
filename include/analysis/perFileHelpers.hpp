#pragma once

#include <regex>
#include <array>
#include <vector>
#include <string>
#include "core/types.hpp"

class TFile;

/// @namespace perFileHelpers
/// @brief Namespace for helper plotting functions producing per-filerelevant statistics
/// for each measurement entry
namespace perFileHelpers {

    struct ColumnsShift {
        int start;
        int end;
        int shift;
    };

    constexpr std::array<ColumnsShift, 1> columnShifts = {{
        // {16, 31, -8} // Shifting strips 16-31 down by 8 to create a continuous range of 0-23 for TB data
    }};

    inline int remapStrip(int rawStrip) {
        for (const auto& col : columnShifts)
            if (rawStrip >= col.start && rawStrip <= col.end)
                return rawStrip + col.shift;
        return rawStrip;
    }

    inline std::pair<double, double> getPercentileRange(
        std::vector<double>& data,
        double lower_p = 0.025,
        double upper_p = 0.975
    ) {
        if (data.empty()) return {-6.5 * TIME_TICK_NS, 6.5 * TIME_TICK_NS};

        std::sort(data.begin(), data.end());
        size_t n = data.size();

        size_t idx_low  = static_cast<size_t>(lower_p * (n - 1));
        size_t idx_high = static_cast<size_t>(upper_p * (n - 1));

        double min_val = data[idx_low];
        double max_val = data[idx_high];

        if (min_val >= max_val) {
            min_val -= TIME_TICK_NS;
            max_val += TIME_TICK_NS;
        }

        return {min_val, max_val};
    }

    inline int extractPrimaryHV(const std::string& filepath) {
        std::regex hv_re("HV_?(\\d+)");
        std::smatch match;
        if (std::regex_search(filepath, match, hv_re)) {
            return std::stoi(match[1].str());
        }
        return -1;
    }

    void plotEfficiencyMap(TFile* input_file, bool use_external_trigger);
    void plotStrip(TFile* input_file);
    void extractBeamSpotRoI(TH2* agg_hist, RoI& roi,
        double beam_threshold = 0.15, double halo_threshold = 0.05);
    void plotDtVsStrip(TFile* input_file, RoI& region_of_interest);
    void plotCS(TFile* input_file);
    void plotToT(TFile* input_file, const RoI& region_of_interest);
    void plotToTVsStrip(TFile* input_file);
    void plotToTVsDtVsStrip(TFile* input_file, const RoI& region_of_interest);
    void plotMultiplicity(TFile* input_file);
    void plotDelay(TFile* input_file);
    void plotToFs(TFile* input_file);
}