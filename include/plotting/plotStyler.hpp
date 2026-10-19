#pragma once

#include <tuple>
#include <string>
#include <vector>

#include "plotTypes.hpp"

class TObject;
class TAxis;
class TCanvas;
class TClass;

namespace PlotterHelpers {
namespace PlotStyler {

    PlotCategory getPlotCategory(const TObject* obj);

    std::optional<std::vector<std::string>> getLabelsFromGroupNames(
        const std::vector<std::string>& group_names);

    std::string matchLabels(const std::string& name, const std::string& pattern);

    struct LVSettings {
        double v_amp_eta1;
        double v_amp_eta2;
        double v_th;
    };
    inline static const std::map<int, LVSettings>& getLVLookupTable() {
        static const std::map<int, LVSettings> lv_table = {
            {1, {1.35, 1.50, 2.20}},
            {2, {1.25, 1.40, 2.20}},
            {3, {1.45, 1.60, 2.20}},
            {4, {1.15, 1.30, 2.20}},
            {5, {1.55, 1.70, 2.20}},
            {6, {1.05, 1.20, 2.20}},
            {7, {1.25, 1.50, 2.20}},
            {78, {1.30, 1.30, 1.80}},
        };
        return lv_table;
    }
    inline static const std::map<double, int>& getFilterLookupTable() {
        static const std::map<double, int> filter_table = {
            {3.3, 538},
            {6.9, 364},
            {10.0, 225},
            {0.0, 0}
        };
        return filter_table;
    }

    std::tuple<std::vector<std::string>, std::string,
    std::string, std::vector<std::string>> compilePlotLabels(
        const std::string& metric_name, TObject* obj);

    void setRange(TObject* obj, TAxis* axis, AxisType axis_type,
                  std::optional<double> default_min,
                  std::optional<double> default_max,
                  const DataCutoffs& cutoffs);

    std::string moveSideContextToAxisTitle(std::string y_axis_label, std::vector<std::string>& title_lines);
    void enforceIntegerMinorTicks(TAxis* axis);
    TGraph* findEfficiencyGraphForObject(TObject* obj);
    void convertMultiGraphToEfficiencyScale(TMultiGraph* mg);
    StylerFnPtr getCustomStyler(PlotCategory category);

    // Styler routines
    void styleEfficiency(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleEfficiencyVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgClusterSizeVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgClusterSizeVsEff(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleRateVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleRateVsEff(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleRateStripVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgToFVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgTRVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgTRVsEff(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleTRStripVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgToTVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgMulVsHV(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleStripDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleStripDistributionCombined(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleRoI(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleCSDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToFDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToFHeatmap(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToTDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToTCorrelation(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToTCombSidesDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToTCombMultiDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleToTCombLayersDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleDelayDistribution(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgMultStrip(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleFracMultStrip(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleAvgDelayStrip(TObject* obj, TCanvas* canvas, TClass* cl);
    void styleDefaultPlot(TObject* obj, TCanvas* canvas, TClass* cl);

} // namespace PlotStyler
} // namespace PlotterHelpers