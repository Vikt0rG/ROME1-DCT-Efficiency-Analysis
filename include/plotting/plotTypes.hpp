#pragma once

#include <optional>
#include <vector>
#include <string>

class TObject;
class TCanvas;
class TClass;

namespace PlotterHelpers {

    enum class PlotCategory {
        Efficiency,
        EfficiencyVsHV,
        StripDistribution,
        StripDistributionCombined,
        ToTDistribution,
        ToTCombinedDistribution,
        ToFDistribution,
        ToFHeatmap,
        AvgToFVsHV,
        TimeResolutionVsHV,
        MeanClusterSizeVsHV,
        RateVsHV,
        AvgToTVsHV,
        AvgMultVsHV,
        Default = -1
    };

    enum class AxisType { X, Y, Z };

    struct DataCutoffs {
        std::optional<double> x_min = std::nullopt;
        std::optional<double> x_max = std::nullopt;
        std::optional<double> y_min = std::nullopt;
        std::optional<double> y_max = std::nullopt;
    };

    using StylerFnPtr = void(*)(TObject*, TCanvas*, TClass*);

} // namespace PlotterHelpers