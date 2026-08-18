#pragma once

#include <tuple>
#include <string>

#include "plotTypes.hpp"

class TObject;
class TAxis;
class TCanvas;
class TClass;

namespace PlotterHelpers {
namespace PlotStyler {

    PlotCategory getPlotCategory(const TObject* obj);

    std::tuple<std::string, std::string, std::string, std::vector<std::string>> compilePlotLabels(
        const std::string& metric_name, TObject* obj);

    void setRange(TObject* obj, TAxis* axis, AxisType axis_type,
                  std::optional<double> default_min,
                  std::optional<double> default_max,
                  const DataCutoffs& cutoffs);

    void enforceIntegerMinorTicks(TAxis* axis);

    StylerFnPtr getCustomStyler(PlotCategory category);

    // Styler routines
    void styleEfficiencyVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgClusterSizeVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgToFVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgToTVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgMulVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleStripDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleStripDistributionCombined(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToFDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToFHeatmap(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToTDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToTCombinedDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleDefaultPlot(TObject* obj, TCanvas* canvas, TClass* cl);

} // namespace PlotStyler
} // namespace PlotterHelpers