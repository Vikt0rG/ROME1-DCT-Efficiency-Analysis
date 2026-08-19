#include <cmath>
#include <utility>
#include <optional>

#include <iostream>
#include <regex>

#include <TObject.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TAxis.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>
#include <TLatex.h>
#include <TColor.h>
#include <TVirtualPad.h>
#include <TMath.h>

#include "plotStyler.hpp"
#include "core/constants.hpp"
#include "ATLASStyler.hpp"
#include "objects.hpp"

namespace PlotterHelpers {
namespace PlotStyler {

    using namespace ATLASStyler;

    static const std::vector<std::tuple<std::string, TClass*, PlotCategory>> category_map = {
        {"h1d_strip_eta",           TH1::Class(),                  PlotCategory::StripDistribution},
        {"h1d_strip_eta",           THStack::Class(),              PlotCategory::StripDistributionCombined},
        {"h1d_tot_eta",             TH1::Class(),                  PlotCategory::ToTDistribution},
        {"h1d_tot",                 THStack::Class(),              PlotCategory::ToTCombinedDistribution},
        {"h1d_tof_layer",           TH1::Class(),                  PlotCategory::ToFDistribution},
        {"h2d_time_of_flight",      TH2::Class(),                  PlotCategory::ToFHeatmap},
        {"avg_time_of_flight",      TMultiGraph::Class(),          PlotCategory::AvgToFVsHV},
        {"time_resolution",         TMultiGraph::Class(),          PlotCategory::TimeResolutionVsHV},
        {"track_eff",               TGraphAsymmErrors::Class(),    PlotCategory::Efficiency},
        {"eff",                     TGraphAsymmErrors::Class(),    PlotCategory::Efficiency},
        {"track_eff",               TMultiGraph::Class(),          PlotCategory::EfficiencyVsHV},
        {"eff",                     TMultiGraph::Class(),          PlotCategory::EfficiencyVsHV},
        {"avg_cluster_size",        TMultiGraph::Class(),          PlotCategory::MeanClusterSizeVsHV},
        {"rate_eta",                TMultiGraph::Class(),          PlotCategory::RateVsHV},
        {"rate_strips_eta",         TMultiGraph::Class(),          PlotCategory::AvgToTVsHV},
        {"avg_tot",                 TMultiGraph::Class(),          PlotCategory::AvgToTVsHV},
        {"avg_multiplicity",        TMultiGraph::Class(),          PlotCategory::AvgMultVsHV}
    };

    static const std::vector<std::pair<PlotCategory, StylerFnPtr>> styler_map = {
        {PlotCategory::EfficiencyVsHV,              &styleEfficiencyVsHV},
        {PlotCategory::MeanClusterSizeVsHV,         &styleAvgClusterSizeVsHV},
        {PlotCategory::RateVsHV,                    &styleAvgClusterSizeVsHV},
        {PlotCategory::ToFDistribution,             &styleToFDistribution},
        {PlotCategory::ToFHeatmap,                  &styleToFHeatmap},
        {PlotCategory::AvgToFVsHV,                  &styleAvgToFVsHV},
        {PlotCategory::TimeResolutionVsHV,          &styleAvgToFVsHV},
        {PlotCategory::AvgToTVsHV,                  &styleAvgToTVsHV},
        {PlotCategory::AvgMultVsHV,                 &styleAvgMulVsHV},
        {PlotCategory::StripDistribution,           &styleStripDistribution},
        {PlotCategory::StripDistributionCombined,   &styleStripDistributionCombined},
        {PlotCategory::ToTDistribution,             &styleToTDistribution},
        {PlotCategory::ToTCombinedDistribution,     &styleToTCombinedDistribution},
        {PlotCategory::Default,                     &styleDefaultPlot}
    };

    PlotCategory getPlotCategory(const TObject* obj) {
        if (!obj) return PlotCategory::Default;
        std::string name = obj->GetName();
        for (const auto& [token, cl, category] : category_map) {
            if (name.find(token) != std::string::npos && obj->InheritsFrom(cl)) {
                return category;
            }
        }
        return PlotCategory::Default;
    }

    StylerFnPtr getCustomStyler(PlotCategory category) {
        for (const auto& [cat, fn_ptr] : styler_map) {
            if (cat == category) return fn_ptr;
        }
        return nullptr;
    }

    std::tuple<std::string, std::string, std::string, std::vector<std::string>> compilePlotLabels(
        const std::string& metric_name,
        TObject* obj)
    {
        std::string out_title, out_xaxis = "High Voltage [V]", out_yaxis;
        std::vector<std::string> title_parts;

        // -------------------------------------------------------------------------
        // Determine Y-Axis based on metric keywords
        if (metric_name.find("eff") != std::string::npos) {
            out_yaxis = "Efficiency";
        } else if (metric_name.find("cluster_size") != std::string::npos) {
            out_yaxis = "#LTCluster Size#GT [Hits]";
        } else if (metric_name.find("rate") != std::string::npos) {
            out_yaxis = "Rate [Hz/cm^{2}]";
        } else if (metric_name.find("tot") != std::string::npos) {
            out_yaxis = "#LTToT#GT [ns]";
        } else if (metric_name.find("multiplicity") != std::string::npos) {
            out_yaxis = "#LTMultiplicity#GT [Hits]";
        } else if (metric_name.find("time_resolution") != std::string::npos) {
            out_yaxis = "Time Resolution [Ticks]";
        } else if (metric_name.find("time_of_flight") != std::string::npos) {
            out_yaxis = metric_name.find("avg_") != std::string::npos ? "#LTToF#GT [Ticks]" : "Time of Flight [Ticks]";
        } else if (metric_name.find("strip") != std::string::npos) {
            out_xaxis = "Strip Number"; out_yaxis = "Hits";
        } else {
            out_xaxis = "X-Axis [a.u.]"; out_yaxis = "Value [a.u.]";
        }

        // -------------------------------------------------------------------------
        // Build Subtitle Context Pieces

        std::smatch match;

        // A. Heatmap Prefix
        if (metric_name.find("h2d_") == 0) {
            title_parts.push_back("Heatmap");
        }

        // B. Track Reconstruction Context
        static const std::regex reco_re("^(track_)?(avg_tot|avg_multiplicity|eff)_");
        if (std::regex_search(metric_name, match, reco_re)) {
            title_parts.push_back(match[1].matched ? "After Track Reco" : "Before Track Reco");
        }

        // C. Layer / Layer Pair Context
        static const std::regex layer_pair_re("layer_(\\d)_(\\d)");
        static const std::regex single_layer_re("layer(\\d+)");

        // 1. Look for the ToF layer pair (e.g., "layer_0_1")
        if (std::regex_search(metric_name, match, layer_pair_re)) {
            std::string new_title = Form("#it{t}_{Layer %s} #minus #it{t}_{Layer %s}",
                                         match[1].str().c_str(),
                                         match[2].str().c_str());
            title_parts.push_back(new_title);
        }

        // 2. Look for a single layer (e.g., "layer0")
        if (std::regex_search(metric_name, match, single_layer_re)) {
            title_parts.push_back("Layer " + match[1].str());
        }

        // D. Side Context
        static const std::regex side_re("eta1|eta2|_or_|_and_");
        if (std::regex_search(metric_name, match, side_re)) {
            std::string m = match.str(0);
            if (m == "eta1")       title_parts.push_back("Side #eta_{1}");
            else if (m == "eta2")  title_parts.push_back("Side #eta_{2}");
            else if (m == "_or_")  title_parts.push_back("OR(#eta_{1}, #eta_{2})");
            else if (m == "_and_") title_parts.push_back("AND(#eta_{1}, #eta_{2})");
        }

        // E. Trigger Context
        static const std::regex trigger_re("external|rpc");
        if (std::regex_search(metric_name, match, trigger_re)) {
            std::string m = match.str(0);
            if (m == "external") title_parts.push_back("External Trigger");
            else if (m == "rpc") title_parts.push_back("RPC Coincidence");
        }

        // Assemble Title (e.g., "After Track Reco: Layers 0 & 1: Side #eta_{1}")
        for (size_t i = 0; i < title_parts.size(); ++i) {
            out_title += title_parts[i];
            if (i < title_parts.size() - 1) out_title += ": ";
        }

        // -------------------------------------------------------------------------
        /// Build legend entries if object is a container (TMultiGraph or THStack)
        auto matchLabels = [](const std::string& name, const std::string& pattern) -> std::string {
            std::regex re(pattern);
            std::smatch match;
            if (std::regex_search(name, match, re)) return match[1].str();
            return "";
        };

        std::vector<std::string> legend_entries;
        if (TMultiGraph* mg = dynamic_cast<TMultiGraph*>(obj)) {
            TIter next(mg->GetListOfGraphs());
            TObject* obj;

            while ((obj = next())) {

                if (auto gr = dynamic_cast<TGraph*>(obj)) {
                    std::string gr_name = gr->GetTitle();
                    std::string layer = "", strip = "";

                    layer = matchLabels(gr_name, "layer(\\d+)");
                    strip = matchLabels(gr_name, "strip(\\d+)");
                    std::string legend_entry;

                    if (!layer.empty()) {
                        legend_entry = "Layer " + layer;
                    }
                    if (!strip.empty()) {
                        legend_entry = "Strip" + strip;
                    }

                    legend_entries.push_back(legend_entry);
                }
            }
        } else if (dynamic_cast<THStack*>(obj)) {
            if (metric_name.find("tot") != std::string::npos) {
                legend_entries.push_back("Side #eta_{1}");
                legend_entries.push_back("Side #eta_{2}");
            } else if (metric_name.find("strip") != std::string::npos) {

                static const std::regex reco_status("before|after|rejected");

                auto words_begin = std::sregex_iterator(metric_name.begin(), metric_name.end(), reco_status);
                auto words_end = std::sregex_iterator();

                for (std::sregex_iterator i = words_begin; i != words_end; ++i) {
                    std::string status = (*i).str();

                    if (status == "after") {
                        legend_entries.push_back("After Track Reco");
                    } else if (status == "before") {
                        legend_entries.push_back("Before Track Reco");
                    } else if (status == "rejected") {
                        legend_entries.push_back("Rejected");
                    }
                }
            }
        }

        return std::make_tuple(out_title, out_xaxis, out_yaxis, legend_entries);
    }

    void setRange(TObject* obj, TAxis* axis, AxisType axis_type,
        std::optional<double> default_min = std::nullopt,
        std::optional<double> default_max = std::nullopt,
        const DataCutoffs& cutoffs = {}) {

        if (!obj || !axis) return;

        double true_min = INT_MAX;
        double true_max = INT_MIN;
        bool found_valid_points = false;

        // Recursive lambda to dynamically parse any ROOT object for min/max limits
        std::function<void(TObject*)> extractBounds = [&](TObject* current_obj) {
            if (!current_obj) return;

            // CONTAINERS 1: TMultiGraph (Unpack and recurse)
            if (auto mg = dynamic_cast<TMultiGraph*>(current_obj)) {
                if (mg->GetListOfGraphs()) {
                    for (TObject* child : *mg->GetListOfGraphs()) {
                        extractBounds(child);
                    }
                }
            }
            // CONTAINERS 2: THStack (Unpack and recurse)
            else if (auto stack = dynamic_cast<THStack*>(current_obj)) {
                if (stack->GetHists()) {
                    for (TObject* child : *stack->GetHists()) {
                        extractBounds(child);
                    }
                }
            }
            // DATA 1: TGraphs (Handles standard, Errors, and AsymmErrors)
            else if (auto gr = dynamic_cast<TGraph*>(current_obj)) {
                auto* gr_err = dynamic_cast<TGraphErrors*>(gr);
                auto* gr_asymm = dynamic_cast<TGraphAsymmErrors*>(gr);

                int n_points = gr->GetN();
                for (int i = 0; i < n_points; ++i) {
                    double x_val = gr->GetX()[i];
                    double y_val = gr->GetY()[i];

                    // Data cutoffs
                    if (cutoffs.x_min.has_value() && x_val < cutoffs.x_min.value()) continue;
                    if (cutoffs.x_max.has_value() && x_val > cutoffs.x_max.value()) continue;
                    if (cutoffs.y_min.has_value() && y_val < cutoffs.y_min.value()) continue;
                    if (cutoffs.y_max.has_value() && y_val > cutoffs.y_max.value()) continue;

                    double val_low = 0.0, val_high = 0.0;
                    if (axis_type == AxisType::X) {
                        val_low = val_high = x_val;
                        if (gr_asymm) {
                            val_low -= gr_asymm->GetErrorXlow(i);
                            val_high += gr_asymm->GetErrorXhigh(i);
                        } else if (gr_err) {
                            val_low -= gr_err->GetErrorX(i);
                            val_high += gr_err->GetErrorX(i);
                        }
                    } else if (axis_type == AxisType::Y) {
                        val_low = val_high = y_val;
                        if (gr_asymm) {
                            val_low -= gr_asymm->GetErrorYlow(i);
                            val_high += gr_asymm->GetErrorYhigh(i);
                        } else if (gr_err) {
                            val_low -= gr_err->GetErrorY(i);
                            val_high += gr_err->GetErrorY(i);
                        }
                    }

                    if (val_low < true_min) true_min = val_low;
                    if (val_high > true_max) true_max = val_high;
                    found_valid_points = true;
                }
            }
            // DATA 2: Histograms (Handles TH1, TH2, TH3 seamlessly)
            else if (auto h = dynamic_cast<TH1*>(current_obj)) {
                int n_bins_x = h->GetNbinsX();
                int n_bins_y = h->GetNbinsY();
                int n_bins_z = h->GetNbinsZ();

                for (int x = 1; x <= n_bins_x; ++x) {
                    double x_center = h->GetXaxis()->GetBinCenter(x);

                    // X-axis cutoffs
                    if (cutoffs.x_min.has_value() && x_center < cutoffs.x_min.value()) continue;
                    if (cutoffs.x_max.has_value() && x_center > cutoffs.x_max.value()) continue;

                    for (int y = 1; y <= n_bins_y; ++y) {
                        for (int z = 1; z <= n_bins_z; ++z) {
                            int global_bin = h->GetBin(x, y, z);
                            double content = h->GetBinContent(global_bin);
                            double error = h->GetBinError(global_bin);

                            // Skip empty bins
                            if (content == 0 && error == 0) continue; 

                            if (cutoffs.y_min.has_value() && content < cutoffs.y_min.value()) continue;
                            if (cutoffs.y_max.has_value() && content > cutoffs.y_max.value()) continue;

                            double val_low = 0.0, val_high = 0.0;
                            if (axis_type == AxisType::X) {
                                val_low = h->GetXaxis()->GetBinLowEdge(x);
                                val_high = h->GetXaxis()->GetBinUpEdge(x);
                            } else if (axis_type == AxisType::Y) {
                                val_low = content - error;
                                val_high = content + error;
                            }

                            if (val_low < true_min) true_min = val_low;
                            if (val_high > true_max) true_max = val_high;
                            found_valid_points = true;
                        }
                    }
                }
            }
        };

        // Start the recursive extraction from the provided object
        extractBounds(obj);

        // Apply dynamic margins if valid points were found
        if (found_valid_points && true_max >= true_min) {
            double safety_buffer = (true_max > true_min) ? (true_max - true_min) * 0.05 : 0.05;

            double dynamic_min = true_min - safety_buffer;
            double dynamic_max = true_max + safety_buffer;

            // Clamping only happens if a default_min is provided
            if (axis_type == AxisType::Y && default_min.has_value() && dynamic_min < default_min.value()) {
                dynamic_min = default_min.value();
            }

            if (axis_type == AxisType::Y && default_max.has_value() && dynamic_max > default_max.value()) {
                dynamic_max = default_max.value();
            }

            axis->SetLimits(dynamic_min, dynamic_max);
            axis->SetRangeUser(dynamic_min, dynamic_max);
        } else {
            if (default_min.has_value() && default_max.has_value()) {
                axis->SetLimits(default_min.value(), default_max.value());
                axis->SetRangeUser(default_min.value(), default_max.value());
            }
        }
    }

    void enforceIntegerMinorTicks(TAxis* axis) {
        if (!axis) return;

        int current_ndiv = axis->GetNdivisions();
        int n1 = current_ndiv % 100;         // Major divisions
        int n2 = (current_ndiv / 100) % 100; // Minor divisions
        int n3 = current_ndiv / 10000;       // Tertiary divisions

        if (n1 == 0) n1 = 1; // Safety against div-by-zero
        double range = axis->GetXmax() - axis->GetXmin();
        double min_major_step = range / n1;

        if (min_major_step / n2 < 1.0) {
            n2 = static_cast<int>(std::floor(min_major_step));
            if (n2 < 1) n2 = 1;
            axis->SetNdivisions(n1 + 100 * n2 + 10000 * n3, kTRUE);
        }
    };

    void styleEfficiencyVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        // Extract title and axis labels from the object's title string
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        obj->Draw("APE0");

        // Set axis ranges and labels
        if (mg && mg->GetHistogram()) {

            // Setup X Axis
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, 0.0, 10e3, {.x_min = 4500.0});
                xAxis->SetTitle(x_label.c_str());
            }

            // Setup Y Axis (Default fallback 0.0 to 1.0, obeying the same X-floor)
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, 0.0, 1.0);
                yAxis->SetTitle(y_label.c_str());
            }
        }

        // Set color and marker style for the graphs in the multigraph
        const std::vector<Color_t> palette = {kAzure + 2, kGreen + 2, kOrange + 10};
        if (mg && mg->GetListOfGraphs()) {
            TIter next(mg->GetListOfGraphs());
            TObject* gr_obj;
            int color_idx = 0;
            while ((gr_obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                    Color_t color = palette[color_idx % palette.size()];

                    gr->SetMarkerStyle(52);
                    gr->SetMarkerSize(1.8);
                    gr->SetMarkerColor(color);
                    gr->SetLineColor(color);
                    gr->SetLineWidth(1);

                    color_idx++;

                    // Extract the fitted sigmoid function from the graph's list of functions
                    TF1* sigmoid = dynamic_cast<TF1*>(gr->GetListOfFunctions()->First());
                    if (!sigmoid) continue;

                    // Sigmoid fit band calculation
                    double p0 = sigmoid->GetParameter(0);
                    double p1 = sigmoid->GetParameter(1);
                    double p2 = sigmoid->GetParameter(2);

                    double ep0 = sigmoid->GetParError(0);
                    double ep1 = sigmoid->GetParError(1);
                    double ep2 = sigmoid->GetParError(2);

                    int n_band_points = 200;
                    double x_min_band = gr->GetXaxis()->GetXmin();
                    double x_max_band = gr->GetXaxis()->GetXmax();
                    double step = (x_max_band - x_min_band) / n_band_points;

                    TGraphErrors* fit_band = new TGraphErrors(n_band_points);
                    double sigma_multiplier = 3.0;

                    for (int i = 0; i < n_band_points; ++i) {
                        double x = x_min_band + i * step;
                        double y = sigmoid->Eval(x);

                        double exponent = -p1 * (x - p2);
                        double dy = 0.0;

                        if (exponent > -50.0 && exponent < 50.0) {
                            double E = TMath::Exp(exponent);
                            double D = 1.0 + E;

                            double df_dp0 = (p0 != 0) ? y / p0 : 0;
                            double df_dp1 = y * (x - p2) * E / D;
                            double df_dp2 = -y * p1 * E / D;

                            dy = sigma_multiplier * std::sqrt(std::pow(df_dp0 * ep0, 2) +
                                                              std::pow(df_dp1 * ep1, 2) +
                                                              std::pow(df_dp2 * ep2, 2));
                        }

                        fit_band->SetPoint(i, x, y);
                        fit_band->SetPointError(i, 0, dy);
                    }

                    fit_band->SetFillColorAlpha(color, 0.3);
                    fit_band->SetLineColor(color);
                    fit_band->SetLineWidth(0);
                    fit_band->Draw("E3 SAME");
                    fit_band->SetBit(kCanDelete);

                    sigmoid->SetLineColor(color);
                    sigmoid->SetLineWidth(2);
                    sigmoid->Draw("SAME");

                    sigmoid->SetBit(kCanDelete);
                }
            }
        }

        if (auto named_obj = dynamic_cast<TNamed*>(obj)) {
            named_obj->SetTitle(title.c_str());
        }

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            plot_title,               // Title string
            12,                       // Alignment
            kWhite, 0.70,             // semi-transparent white background
            kBlack, 1,                // Black 1px border line
            0.01                      // Inner padding
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
        drawATLASLegend(obj, legend_entries, 0.18, legend_y, 13);

        canvas->Modified();
        canvas->Update();
    }

    void styleAvgClusterSizeVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        // Extract title and axis labels from the object's title string
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        obj->Draw("AP0Z");

        // Set axis ranges and labels
        if (mg && mg->GetHistogram()) {

            // Setup X Axis
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 4500.0});
                xAxis->SetTitle(x_label.c_str());
            }

            // Setup Y Axis (Default fallback 0.0 to 1.0, obeying the same X-floor)
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y);
                yAxis->SetTitle(y_label.c_str());
            }
        }

        // Set color and marker style for the graphs in the multigraph
        const std::vector<Color_t> palette = {kAzure + 2, kGreen + 2, kOrange + 10};
        if (mg && mg->GetListOfGraphs()) {
            TIter next(mg->GetListOfGraphs());
            TObject* gr_obj;
            int color_idx = 0;
            while ((gr_obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                    Color_t color = palette[color_idx % palette.size()];

                    gr->SetMarkerStyle(52);
                    gr->SetMarkerSize(1.8);
                    gr->SetMarkerColor(color);
                    gr->SetLineColor(color);
                    gr->SetLineWidth(1);

                    color_idx++;
                }
            }
        }

        if (auto named_obj = dynamic_cast<TNamed*>(obj)) {
            named_obj->SetTitle(title.c_str());
        }

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            plot_title,               // Title string
            12,                       // Alignment
            kWhite, 0.70,             // semi-transparent white background
            kBlack, 1,                // Black 1px border line
            0.01                      // Inner padding
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
        drawATLASLegend(obj, legend_entries, 0.18, legend_y, 13);

        canvas->Modified();
        canvas->Update();
    }

    void styleAvgToFVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {
        // Extract title and axis labels from the object's title string
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        obj->Draw("AP0Z");

        // Set axis ranges and labels
        if (mg && mg->GetHistogram()) {

            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 5200.0});
                xAxis->SetTitle(x_label.c_str());
            }

            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 5200.0});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        // Set color and marker style for the graphs in the multigraph
        const std::vector<Color_t> palette = {kAzure + 2, kGreen + 2, kOrange + 10};
        if (mg && mg->GetListOfGraphs()) {
            TIter next(mg->GetListOfGraphs());
            TObject* gr_obj;
            int color_idx = 0;
            while ((gr_obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                    Color_t color = palette[color_idx % palette.size()];

                    gr->SetMarkerStyle(52);
                    gr->SetMarkerSize(1.8);
                    gr->SetMarkerColor(color);
                    gr->SetLineColor(color);
                    gr->SetLineWidth(1);

                    color_idx++;
                }
            }
        }

        if (auto named_obj = dynamic_cast<TNamed*>(obj)) {
            named_obj->SetTitle(title.c_str());
        }

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.03,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            plot_title,               // Title string
            32,                       // Alignment
            kWhite, 0.70,             // semi-transparent white background
            kBlack, 1,                // Black 1px border line
            0.01                      // Inner padding
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header->GetY1NDC() - 0.02;
        int alignment = 33;
        drawATLASLegend(obj, legend_entries, ndc_x0, legend_y, alignment);

        canvas->Modified();
        canvas->Update();
    }

    void styleAvgToTVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        // Differentiate between a 24-strip plot and a 3-layer plot
        bool is_strip_plot = (mg && mg->GetListOfGraphs() && mg->GetListOfGraphs()->GetSize() > 10);

        if (is_strip_plot) gStyle->SetPalette(kViridis);

        obj->Draw("APZ");

        if (mg && mg->GetHistogram()) {
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 4500.0});
                xAxis->SetTitle(x_label.c_str());
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt);
                yAxis->SetTitle(y_label.c_str());
            }
        }

        if (auto named_obj = dynamic_cast<TNamed*>(obj)) {
            named_obj->SetTitle(title.c_str());
        }

        applyATLASStyle(obj, canvas);

        // Add a strip colorbar
        if (is_strip_plot) {
            canvas->SetRightMargin(0.16);

            int n_colors = TColor::GetNumberOfColors();
            int max_strip = STRIPS_PER_LAYER - 1;

            // Recolor the 24 graphs to match the continuous palette
            TIter next(mg->GetListOfGraphs());
            TObject* gr_obj;
            while ((gr_obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                    int strip_idx = 0;
                    std::smatch match;
                    std::string gr_title = gr->GetTitle();

                    // Extract the strip number from the title (e.g. "Strip 5")
                    if (std::regex_search(gr_title, match, std::regex("Strip (\\d+)"))) {
                        strip_idx = std::stoi(match[1].str());
                    }
                    strip_idx = std::max(0, std::min(strip_idx, max_strip)); // Safety clamp

                    // Map the strip number [0, 23] to the palette index [0, 255]
                    int color_idx = TColor::GetColorPalette((strip_idx * (n_colors - 1)) / max_strip);

                    gr->SetMarkerColor(color_idx);
                    gr->SetMarkerStyle(70);
                    gr->SetMarkerSize(1.8);
                    gr->SetLineColor(color_idx);
                    gr->SetLineWidth(2.0);
                }
            }

            // Create a dummy histogram specifically to draw the Z-axis (Colorbar)
            TH2D* dummy_z = new TH2D(Form("dummy_z_%p", mg), "", 1, -2000, -1000, 1, -2000, -1000);
            dummy_z->SetDirectory(nullptr);
            dummy_z->SetBinContent(1, 1, 0.0);
            dummy_z->SetMinimum(0);
            dummy_z->SetMaximum(STRIPS_PER_LAYER);
            dummy_z->SetContour(STRIPS_PER_LAYER);

            TAxis* zAxis = dummy_z->GetZaxis();
            zAxis->SetTitle("Strip Number");
            zAxis->SetTitleOffset(1.0);
            zAxis->SetTitleSize(0.05);
            zAxis->SetLabelSize(0.04);
            zAxis->SetNdivisions(6, 4, 0, kFALSE);

            dummy_z->Draw("COL Z SAME");
        }

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            plot_title,
            12,
            kWhite, 0.70,
            kBlack, 1,
            0.01
        );

        canvas->Modified();
        canvas->Update();

        if (!is_strip_plot) {
            double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
            TLegend* leg = drawATLASLegend(obj, legend_entries, 0.18, legend_y, 13);
            if (leg) {
                leg->SetBorderSize(1);
                leg->SetLineWidth(1);
                leg->SetLineColor(kBlack);
            }
        }

        canvas->Modified();
        canvas->Update();
    }

    void styleAvgMulVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        // Differentiate between a 24-strip plot and a 3-layer plot
        bool is_strip_plot = (mg && mg->GetListOfGraphs() && mg->GetListOfGraphs()->GetSize() > 10);

        if (is_strip_plot) gStyle->SetPalette(kViridis);

        obj->Draw("APZ");

        if (mg && mg->GetHistogram()) {
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 4500.0});
                xAxis->SetTitle(x_label.c_str());
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 4500.0});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        if (auto named_obj = dynamic_cast<TNamed*>(obj)) {
            named_obj->SetTitle(title.c_str());
        }

        applyATLASStyle(obj, canvas);

        // Add a strip colorbar
        if (is_strip_plot) {
            canvas->SetRightMargin(0.16);

            int n_colors = TColor::GetNumberOfColors();
            int max_strip = STRIPS_PER_LAYER - 1;

            // Recolor the 24 graphs to match the continuous palette
            TIter next(mg->GetListOfGraphs());
            TObject* gr_obj;
            while ((gr_obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                    int strip_idx = 0;
                    std::smatch match;
                    std::string gr_title = gr->GetTitle();

                    // Extract the strip number from the title (e.g. "Strip 5")
                    if (std::regex_search(gr_title, match, std::regex("Strip (\\d+)"))) {
                        strip_idx = std::stoi(match[1].str());
                    }
                    strip_idx = std::max(0, std::min(strip_idx, max_strip)); // Safety clamp

                    // Map the strip number [0, 23] to the palette index [0, 255]
                    int color_idx = TColor::GetColorPalette((strip_idx * (n_colors - 1)) / max_strip);

                    gr->SetMarkerColor(color_idx);
                    gr->SetMarkerStyle(70);
                    gr->SetMarkerSize(1.8);
                    gr->SetLineColor(color_idx);
                    gr->SetLineWidth(2.0);
                }
            }

            // Create a dummy histogram specifically to draw the Z-axis (Colorbar)
            TH2D* dummy_z = new TH2D(Form("dummy_z_%p", mg), "", 1, -2000, -1000, 1, -2000, -1000);
            dummy_z->SetDirectory(nullptr);
            dummy_z->SetBinContent(1, 1, 0.0);
            dummy_z->SetMinimum(0);
            dummy_z->SetMaximum(STRIPS_PER_LAYER);
            dummy_z->SetContour(STRIPS_PER_LAYER);

            TAxis* zAxis = dummy_z->GetZaxis();
            zAxis->SetTitle("Strip Number");
            zAxis->SetTitleOffset(1.0);
            zAxis->SetTitleSize(0.05);
            zAxis->SetLabelSize(0.04);
            zAxis->SetNdivisions(6, 4, 0, kFALSE);

            dummy_z->Draw("COL Z SAME");
        }

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            plot_title,
            12,
            kWhite, 0.70,
            kBlack, 1,
            0.01
        );

        canvas->Modified();
        canvas->Update();

        if (!is_strip_plot) {
            double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
            TLegend* leg = drawATLASLegend(obj, legend_entries, 0.18, legend_y, 13);
            if (leg) {
                leg->SetBorderSize(1);
                leg->SetLineWidth(1);
                leg->SetLineColor(kBlack);
            }
        }

        canvas->Modified();
        canvas->Update();
    }

    void styleStripDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {

        auto h1 = dynamic_cast<TH1*>(obj);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(2.0);
        h1->SetLineStyle(1);

        // Format: TColor::GetColorTransparent(Color_Index, Alpha_Opacity_From_0_to_1)
        Int_t light_blue_transparent = TColor::GetColorTransparent(kAzure + 7, 0.30);
        h1->SetFillColor(light_blue_transparent);
        h1->SetFillStyle(1001); // 1001 = Solid fill style

        h1->Draw("HIST");

        applyATLASStyle(obj, canvas);

        drawATLASLabel(0.21, 0.86, "Work in Progress");
        drawPlotTitle(obj, 0.21, 0.82);
    }

    void styleStripDistributionCombined(TObject* obj, TCanvas* canvas, TClass* cl) {

        std::cout << "Styling THStack object: " << obj->GetName() << std::endl;
        auto stack = dynamic_cast<THStack*>(obj);
        if (!stack) return;

        stack->Draw("nostack hist");

        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetName(), stack);
        if (auto named_obj = dynamic_cast<TNamed*>(obj)) {
            named_obj->SetTitle(title.c_str());
        }

        TIter next(stack->GetHists());
        TH1* hist = nullptr;
        int index = 0;

        while ((hist = static_cast<TH1*>(next()))) {
            Color_t base_color = (index == 0) ? kAzure - 3 : kRed + 1;

            hist->SetLineColor(base_color);
            hist->SetLineWidth(2);
            hist->SetLineStyle(1);

            Int_t trans_color = TColor::GetColorTransparent(base_color, 0.30);
            hist->SetFillColor(trans_color);
            hist->SetFillStyle(1001);

            index++;
        }

        if (stack->GetXaxis()) stack->GetXaxis()->SetTitle(x_label.c_str());
        if (stack->GetYaxis()) stack->GetYaxis()->SetTitle(y_label.c_str());
        applyATLASStyle(obj, canvas);

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.09,
            "Work in Progress",
            plot_title,
            12,
            kWhite, 0.70,
            kBlack, 0,
            0.01
        );

        drawATLASLegend(obj, legend_entries, ndc_x0 + 0.03, ndc_y0 - 0.2, 13);

        canvas->Modified();
        canvas->Update();
    }

    void styleToFDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto h1 = dynamic_cast<TH1*>(obj);
        if (!h1) return;

        // Force a larger x and y range to accommodate the legend and confidence band
        double x_max = h1->GetXaxis()->GetXmax();
        double y_max = h1->GetMaximum() * 1.07;

        TH1F* frame = canvas->DrawFrame(-11.0, 0.0, x_max, y_max);
        frame->SetTitle(h1->GetTitle());
        frame->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        frame->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());

        applyATLASStyle(frame, canvas);

        h1->SetFillColor(kOrange - 2);
        h1->SetFillStyle(1001);
        h1->SetLineColor(kOrange + 7);
        h1->SetLineWidth(2);
        h1->SetLineStyle(1);

        h1->Draw("HIST SAME");

        if (TF1* fit = h1->GetFunction("gaus")) {
            double p0 = fit->GetParameter(0);
            double p1 = fit->GetParameter(1);
            double p2 = fit->GetParameter(2);

            double ep0 = fit->GetParError(0);
            double ep1 = fit->GetParError(1);
            double ep2 = fit->GetParError(2);

            int n_points = 200;
            double x_min = h1->GetXaxis()->GetXmin();
            double step = (x_max - x_min) / n_points;

            TGraphErrors* fit_band = new TGraphErrors(n_points);

            double sigma_multiplier = 3.0;
            for (int i = 0; i < n_points; ++i) {
                double x = x_min + i * step;
                double y = fit->Eval(x);

                double df_dp0 = (p0 != 0) ? y / p0 : 0;
                double df_dp1 = (p2 != 0) ? y * (x - p1) / (p2 * p2) : 0;
                double df_dp2 = (p2 != 0) ? y * ((x - p1) * (x - p1)) / (p2 * p2 * p2) : 0;

                double dy = sigma_multiplier * std::sqrt(std::pow(df_dp0 * ep0, 2) +
                                            std::pow(df_dp1 * ep1, 2) +
                                            std::pow(df_dp2 * ep2, 2));

                fit_band->SetPoint(i, x, y);
                fit_band->SetPointError(i, 0, dy);
            }

            fit_band->SetFillColorAlpha(kAzure - 2, 0.4);
            fit_band->SetLineColor(kBlue + 3);
            fit_band->SetLineWidth(2);
            fit_band->SetLineStyle(1);

            fit_band->Draw("E3 SAME");

            fit->SetLineColor(kBlue + 3);
            fit->SetLineWidth(2);
            fit->SetLineStyle(1);
            fit->Draw("SAME");

            fit_band->SetBit(kCanDelete);

            double ndc_x0 = canvas->GetLeftMargin();
            double ndc_y0 = 1.0 - canvas->GetTopMargin();

            TLegend* fit_legend = new TLegend(ndc_x0 + 0.03, ndc_y0 - 0.30, ndc_x0 + 0.35, ndc_y0 - 0.15);
            fit_legend->SetTextAlign(12);
            fit_legend->SetBorderSize(0);
            fit_legend->SetFillStyle(0);
            fit_legend->SetTextFont(42);
            fit_legend->SetTextSize(0.035);

            fit_legend->AddEntry(fit_band, "Gaussian Fit #pm 3#sigma Conf.", "fl");
            fit_legend->AddEntry((TObject*)nullptr, Form("#mu = %.2f #pm %.2f", p1, ep1), "");
            fit_legend->AddEntry((TObject*)nullptr, Form("#sigma = %.2f #pm %.2f", p2, ep2), "");

            fit_legend->Draw();

            std::string plot_title = obj ? obj->GetTitle() : "";
            drawATLASHeaderBlock(
                ndc_x0 + 0.03, ndc_y0 - 0.09,
                "Work in Progress",
                plot_title,
                12,
                kWhite, 0.70,
                kBlack, 0,
                0.01
            );
        }
    }

    void styleToFHeatmap(TObject* obj, TCanvas* canvas, TClass* cl) {

        auto h2 = dynamic_cast<TH2*>(obj);
        h2->Draw("COLZ");

        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), h2);
        if (auto named_obj = dynamic_cast<TNamed*>(obj)) {
            named_obj->SetTitle(title.c_str());
        }

        applyATLASStyle(obj, canvas);

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            plot_title,
            12,
            kWhite, 0.70,
            kBlack, 0,
            0.01
        );
    }

    void styleToTDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {

        auto h1 = dynamic_cast<TH1*>(obj);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(2.0);
        h1->SetLineStyle(1);

        Int_t light_green_transparent = TColor::GetColorTransparent(kGreen + 2, 0.30);
        h1->SetFillColor(light_green_transparent);
        h1->SetFillStyle(1001);

        h1->Draw("HIST");

        applyATLASStyle(obj, canvas);

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            plot_title,
            32,
            kWhite, 0.00,
            kBlack, 0,
            0.01
        );
    }

    void styleToTCombinedDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);
        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), stack);
        if (!stack) return;

        TList* hist_list = stack->GetHists();
        if (!hist_list || hist_list->GetSize() < 2) return; // We need at least 2 histograms to compute a ratio

        canvas->SetWindowSize(800, 840);
        canvas->SetCanvasSize(800, 840);

        // Create the Top and Bottom Pads (Ratio Layout)
        canvas->cd();
        
        // Top Pad: Occupies upper 70% of the canvas
        TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
        pad1->SetLeftMargin(0.16);
        pad1->SetRightMargin(0.05);
        pad1->SetTopMargin(0.07);
        pad1->SetBottomMargin(0.03);
        pad1->Draw();

        // Bottom Pad: Occupies lower 30% of the canvas
        TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.30);
        pad2->SetLeftMargin(0.16);
        pad2->SetRightMargin(0.05);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.35);
        pad2->SetTicks(1, 1);
        pad2->Draw();

        // Style & Extract Histograms from Stack
        std::vector<std::pair<double, Color_t>> mean_line_data;
        TH1* h1 = nullptr;
        TH1* h2 = nullptr;

        TIter next(hist_list);
        TH1* hist = nullptr;
        int index = 0;

        while ((hist = static_cast<TH1*>(next()))) {
            Color_t base_color = (index == 0) ? kBlue - 2 : kRed - 3;
            mean_line_data.push_back({hist->GetMean(), base_color});

            // Store references to the first two histograms for our ratio
            if (index == 0) h1 = hist;
            if (index == 1) h2 = hist;

            // Apply standard styling
            hist->SetLineColor(base_color);
            hist->SetLineWidth(2);
            hist->SetLineStyle(1);

            Int_t trans_color = TColor::GetColorTransparent(base_color, 0.30);
            hist->SetFillColor(trans_color);
            hist->SetFillStyle(1001);

            index++;
        }

        // Draw Top Pad (Main Stack Distributions)
        pad1->cd();
        applyATLASStyle(obj, pad1);
        stack->Draw("nostack hist");

        // Remove X-axis labels/titles from the top plot to avoid overlaps
        if (stack->GetXaxis()) {
            stack->GetXaxis()->SetLabelSize(0);
            stack->GetXaxis()->SetTitleSize(0);
        }
        if (stack->GetYaxis()) {
            stack->GetYaxis()->SetTitle("Hits");
            stack->GetYaxis()->SetTitleSize(0.06);
            stack->GetYaxis()->SetLabelSize(0.05);
            stack->GetYaxis()->SetTitleOffset(1.2);
        }

        pad1->Modified();
        pad1->Update();

        double x_max = pad1->GetUxmax();
        double y_min = pad1->GetUymin();
        double y_max = pad1->GetUymax();

        Objects::line(pad1, TOT_ROI_MAX, TOT_ROI_MAX, y_min, y_max, kBlack, 9, 1);
        Objects::hatchedRegion(pad1, TOT_ROI_MAX, x_max, y_min, y_max, 3244);

        // Draw Vertical Mean Lines on Top Pad
        for (const auto& data : mean_line_data) {
            double mean_x = data.first;
            Color_t color = data.second;

            TLine* mean_line = new TLine(mean_x, y_min, mean_x, y_max);
            mean_line->SetLineColor(color);
            mean_line->SetLineWidth(2);
            mean_line->SetLineStyle(11);
            mean_line->Draw();
        }

        double ndc_x0 = 1.0 - pad1->GetRightMargin();
        double ndc_y0 = 1.0 - pad1->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.05,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            plot_title,               // Title string
            32,                       // Alignment
            kWhite, 0.70f,            // 85% semi-transparent white background
            kBlack, 1,                // Black 1px border line
            0.01                      // Inner padding
        );

        pad1->Modified();
        pad1->Update();

        double legend_y = header ? header->GetY1NDC() - 0.15 : 0.70;
        drawATLASLegend(stack, legend_entries, 0.92, legend_y, 31);

        // Calculate & Draw Bottom Pad (Ratio Plot)
        pad2->cd();

        // Clone h1 to preserve properties, then divide by h2
        TH1* h_ratio = static_cast<TH1*>(h1->Clone("h_ratio"));
        h_ratio->SetDirectory(nullptr);
        h_ratio->Reset();
        h_ratio->Divide(h1, h2, 1.0, 1.0, "B"); // "B" uses binomial errors

        // Stylize Ratio Histogram
        h_ratio->SetLineColor(kBlack);
        h_ratio->SetMarkerColor(kBlack);
        h_ratio->SetLineWidth(2);
        h_ratio->SetMarkerStyle(20);
        h_ratio->SetMarkerSize(0.8);
        h_ratio->SetFillStyle(0);

        h_ratio->Draw("ep"); // Draw with error bars and points

        TAxis* rx = h_ratio->GetXaxis();
        TAxis* ry = h_ratio->GetYaxis();

        if (rx) {
            rx->SetTitle("ToT [ns]");
            rx->SetTitleFont(42);
            rx->SetTitleSize(0.14);  // Scaled up for small pad height
            rx->SetTitleOffset(1.0);

            rx->SetLabelFont(42);
            rx->SetLabelSize(0.11);

            rx->SetNdivisions(505);
            rx->SetTickLength(0.06);
        }

        if (ry) {
            ry->SetTitle("#eta1 / #eta2");
            ry->SetTitleFont(42);
            ry->SetLabelFont(42);
            ry->SetTitleSize(0.14);
            ry->SetLabelSize(0.11);
            ry->SetTitleOffset(0.50);

            ry->SetNdivisions(505);

            h_ratio->SetMinimum(0.0); 
            h_ratio->SetMaximum(2.0);
        }

        // Draw a dashed reference line at Y = 1.0
        double x_min_pad2 = rx ? rx->GetXmin() : 0;
        double x_max_pad2 = rx ? rx->GetXmax() : 100;
        Objects::line(pad2, x_min_pad2, x_max_pad2, 1.0, 1.0, kGray+2, 2, 1);

        pad2->Modified();
        pad2->Update();

        canvas->cd();
        canvas->Modified();
        canvas->Update();
    }

    void styleDefaultPlot(TObject* obj, TCanvas* canvas, TClass* cl) {

        if (cl->InheritsFrom(TH2::Class())) {
            obj->Draw("COLZ");
        } else if (cl->InheritsFrom(TMultiGraph::Class()) || cl->InheritsFrom(TGraphAsymmErrors::Class())) {
            obj->Draw("AP");
            gStyle->SetEndErrorSize(8);
        } else if (cl->InheritsFrom(TGraphErrors::Class())) {
            obj->Draw("APZ");
        } else if (cl->InheritsFrom(TGraph::Class())) {
            obj->Draw("AP");
        } else {
            obj->Draw("E1 X0");
        }

        applyATLASStyle(obj, canvas);

        if (auto h2 = dynamic_cast<TH2*>(obj)) {
            enforceIntegerMinorTicks(h2->GetXaxis());
        } else if (auto h1 = dynamic_cast<TH1*>(obj)) {
            enforceIntegerMinorTicks(h1->GetXaxis());
        }

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            plot_title,
            12,
            kWhite, 0.70,
            kBlack, 1,
            0.01
        );
    }

} // namespace PlotStyler
} // namespace PlotterHelpers