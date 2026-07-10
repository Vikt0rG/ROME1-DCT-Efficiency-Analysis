#pragma once

#include <string>
#include <vector>
#include <tuple>

class TObject;
class TCanvas;
class TClass;
class TPad;
class TH2;
class TMultiGraph;

// ------------------------------------------------------------------------------------------
/// @enum PlotCategory
/// @brief Enum to represent different categories of plots that can be generated
enum class PlotCategory {
    Efficiency,
    EfficiencyVsHV,
    StripDistribution,
    DtVsStrip,
    ToT,
    MeanClusterSizeVsHV,
    NoiseRateVsHV,
    Default
};

// ------------------------------------------------------------------------------------------
namespace PlotterHelpers {
namespace PlotStyler {

    /// @brief Type alias for a function pointer that applies styling to a ROOT object
    using StylerFnPtr = void(*)(TObject*, TCanvas*, TClass*);

    /// @brief Helper function to determine the plot category based on a ROOT object
    /// @param obj A pointer to the ROOT object (e.g., TGraph, TH1, TMultiGraph) to be
    /// categorized
    /// @return The corresponding PlotCategory enum value
    PlotCategory getPlotCategory(const TObject* obj);

    /// @brief Helper function to compile plot labels from a given title string
    /// @param metric_name Metric name string from which to derive plot labels
    /// @return A tuple containing the plot title, x-axis label, and y-axis label
    std::tuple<std::string, std::string, std::string> compilePlotLabels(
        const std::string& metric_name
    );

    /// @brief Helper function to retrieve a custom styler function based on the plot category
    /// @param category PlotCategory enum value for which to retrieve the styler function
    /// @return Function pointer to the corresponding styler function for a given plot category
    StylerFnPtr getCustomStyler(PlotCategory category);

    // --------------------------------------------------------------------------------------
    /// @namespace ATLASStyler
    /// @brief Collection of helper functions for applying ATLAS styling to plots
    namespace ATLASStyler {

        /// @brief Helper function to draw the ATLAS label on a plot
        /// @param x_offset The horizontal offset for the label position
        /// @param y_offset The vertical offset for the label position
        /// @param status The status text to display (e.g., "Work in Progress")
        void drawATLASLabel(double x_offset, double y_offset, const std::string& status);

        /// @brief Helper function to draw a title on a plot
        /// @param obj The ROOT object (e.g., TGraph, TH1) to which the title will be applied
        /// @param x_offset The horizontal offset for the title position
        /// @param y_offset The vertical offset for the title position
        /// @return A boolean indicating whether the title was successfully drawn
        bool drawPlotTitle(TObject* obj, double x_offset, double y_offset);

        /// @brief Helper function to draw a legend on a plot
        /// @param mg A pointer to the TMultiGraph object for which the legend will be drawn
        /// @param x_offset The horizontal offset for the legend position
        /// @param y_offset The vertical offset for the legend position
        void drawATLASLegend(TMultiGraph* mg, double x_offset, double y_offset);

        /// @brief Helper function to adjust the color bar of a 2D histogram dynamically
        /// based on the maximum value in the histogram
        /// @param h2 A pointer to the TH2 histogram object
        /// @param pad A pointer to the TPad object on which the histogram is drawn
        void adjustDynamicCB(TH2* h2, TPad* pad);

        /// @brief Helper function to apply ATLAS styling to a given plot object and canvas
        /// @param obj The ROOT object (e.g., TGraph, TH1) to which the style will be applied
        /// @param pad A pointer to the TPad object on which the object is drawn (optional)
        void applyATLASStyle(TObject* obj, TPad* pad = nullptr);
    }   // namespace ATLASStyler
    // --------------------------------------------------------------------------------------

    /// @brief Helper function to style efficiency vs HV plots
    /// @param obj The ROOT object (e.g., TGraph, TH1) to which the style will be applied
    /// @param canvas The TCanvas on which the object is drawn
    /// @param class_type The TClass pointer representing the type of the ROOT object
    void styleEfficiencyVsHV(TObject* obj, TCanvas* canvas, TClass* class_type);

    /// @brief Helper function to style mean cluster size vs HV plots
    /// @param obj The ROOT object (e.g., TGraph, TH1) to which the style will be applied
    /// @param canvas The TCanvas on which the object is drawn
    /// @param class_type The TClass pointer representing the type of the ROOT object
    void styleMeanClusterSizeVsHV(TObject* obj, TCanvas* canvas, TClass* class_type);

    /// @brief Helper function to style strip distribution plots
    /// @param obj The ROOT object (e.g., TGraph, TH1) to which the style will be applied
    /// @param canvas The TCanvas on which the object is drawn
    /// @param class_type The TClass pointer representing the type of the ROOT object
    void styleStripDistribution(TObject* obj, TCanvas* canvas, TClass* class_type);

    /// @brief Default styling function for generic plots
    /// @param obj The ROOT object (e.g., TGraph, TH1) to which the style will be applied
    /// @param canvas The TCanvas on which the object is drawn
    /// @param class_type The TClass pointer representing the type of the ROOT object
    void styleDefaultPlot(TObject* obj, TCanvas* canvas, TClass* class_type);

}   // namespace PlotStyler
}   // namespace PlotterHelpers