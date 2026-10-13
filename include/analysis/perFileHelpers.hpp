#pragma once

#include <regex>
#include <array>
#include <vector>
#include <string>
#include "core/types.hpp"

class TFile;
class TH2;

/// @brief Validate that all provided TTrees have the same number of entries
/// @param named_trees A list of pairs containing tree names and their corresponding TTree pointers
/// @throws std::runtime_error if the number of entries in any tree differs from the reference tree
inline void requireEqualTreeEntries(std::initializer_list<std::pair<std::string, TTree*>> named_trees) {
    if (named_trees.size() <= 1) return;

    auto it = named_trees.begin();
    std::string ref_name = it->first;
    TTree* ref_tree = it->second;

    if (!ref_tree || ref_tree->IsZombie()) {
        throw std::runtime_error("Invalid or null TTree provided for '" + ref_name + "'.");
    }

    Long64_t ref_entries = ref_tree->GetEntries();

    for (++it; it != named_trees.end(); ++it) {
        std::string tree_name = it->first;
        TTree* tree = it->second;

        if (!tree || tree->IsZombie()) {
            throw std::runtime_error("Invalid or null TTree provided for '" + tree_name + "'.");
        }

        Long64_t entries = tree->GetEntries();
        if (entries != ref_entries) {
            throw std::runtime_error(
                "Tree entry mismatch detected! '" + ref_name + "' has " +
                std::to_string(ref_entries) + " entries, but '" + tree_name +
                "' has " + std::to_string(entries) + " entries."
            );
        }
    }
}

/// @brief Validate that all provided data containers have the same size
/// @param named_sizes A list of pairs containing data container names and their corresponding sizes
inline void requireEqualDataSizes(std::initializer_list<std::pair<std::string, size_t>> named_sizes) {
    if (named_sizes.size() <= 1) return;

    auto it = named_sizes.begin();
    const std::string& expected_name = it->first;
    size_t expected_size = it->second;

    for (++it; it != named_sizes.end(); ++it) {
        if (it->second != expected_size) {
            throw std::runtime_error(
                "Data size mismatch detected! '" + expected_name +
                "' has size " + std::to_string(expected_size) +
                ", but '" + it->first +
                "' has size " + std::to_string(it->second) + "."
            );
        }
    }
}

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
        {16, 31, -8} // Shifting strips 16-31 down by 8 to create a continuous range of 0-23 for TB data
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