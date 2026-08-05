#include <cmath>
#include <utility>

#include <iostream>
#include <regex>

#include <TObject.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TAxis.h>
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

#include "plotting/plotStyler.hpp"
#include "core/constants.hpp"

namespace PlotterHelpers {
namespace PlotStyler {

    using namespace ATLASStyler;

    /// @brief Data-driven mapping of metric name prefixes to corresponding plot categories
    static const std::vector<std::tuple<std::string, TClass*, PlotCategory>> category_map = {
        {"h1d_strip_eta",           TH1::Class(),                  PlotCategory::StripDistribution},
        {"h1d_tot_eta",             TH1::Class(),                  PlotCategory::ToTDistribution},
        {"h1d_tot",                 THStack::Class(),              PlotCategory::ToTCombinedDistribution},
        {"track_eff",               TGraphAsymmErrors::Class(),    PlotCategory::Efficiency},
        {"eff",                     TGraphAsymmErrors::Class(),    PlotCategory::Efficiency},
        {"track_eff",               TMultiGraph::Class(),          PlotCategory::EfficiencyVsHV},
        {"eff",                     TMultiGraph::Class(),          PlotCategory::EfficiencyVsHV},
        {"avg_cluster_size",        TMultiGraph::Class(),          PlotCategory::MeanClusterSizeVsHV},
        {"noise_rate_eta",          TMultiGraph::Class(),          PlotCategory::NoiseRateVsHV},
        {"noise_rate_strips_eta",   TMultiGraph::Class(),          PlotCategory::AvgToTVsHV},
        {"avg_tot",                 TMultiGraph::Class(),          PlotCategory::AvgToTVsHV},
        {"avg_multiplicity",        TMultiGraph::Class(),          PlotCategory::AvgMultVsHV}
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

    std::tuple<std::string, std::string, std::string, std::vector<std::string>> compilePlotLabels(
        const std::string& metric_name,
        TObject* obj)
    {
        std::string out_title = "", out_xaxis = "", out_yaxis = "";
        std::string side = "", trigger = "";

        // Match metric name patterns to determine axis labels and (sub)titles
        if (metric_name.rfind("track_eff_", 0) == 0 || metric_name.rfind("eff_", 0) == 0) {
            out_xaxis = "High Voltage [V]"; out_yaxis = "Efficiency";
        } else if (metric_name.rfind("avg_cluster_size_", 0) == 0) {
            out_xaxis = "High Voltage [V]"; out_yaxis = "#LTCluster Size#GT [Hits]";
        } else if (metric_name.rfind("noise_rate_", 0) == 0) {
            out_xaxis = "High Voltage [V]"; out_yaxis = "Noise Rate [Hz/cm^{2}]";
        } else if (metric_name.rfind("track_avg_tot_", 0) == 0 || metric_name.rfind("avg_tot_", 0) == 0) {
            out_xaxis = "High Voltage [V]"; out_yaxis = "#LTToT#GT [ns]";
        } else if (metric_name.rfind( "track_avg_multiplicity_", 0) == 0 || metric_name.rfind( "avg_multiplicity_", 0) == 0) {
            out_xaxis = "High Voltage [V]"; out_yaxis = "#LTMultiplicity#GT [Hits]";
        }

        auto matchLabels = [](const std::string& name, const std::string& pattern) -> std::string {
            std::regex re(pattern);
            std::smatch match;
            if (std::regex_search(name, match, re)) {
                return match[1].str();
            }
            return "";
        };

        /// Build legend entries if object is a container (TMultiGraph or THStack)
        std::vector<std::string> legend_entries;
        if (TMultiGraph* mg = dynamic_cast<TMultiGraph*>(obj)) {
            TIter next(mg->GetListOfGraphs());
            TObject* obj;
            while ((obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(obj)) {
                    std::string gr_name = gr->GetName();
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
            legend_entries.push_back("Side #eta_{1}");
            legend_entries.push_back("Side #eta_{2}");
        }

        // Get side and trigger information directly from metric name
        if (metric_name.find("eta1") != std::string::npos)       side = "#eta_{1} Side";
        else if (metric_name.find("eta2") != std::string::npos)  side = "#eta_{2} Side";
        else if (metric_name.find("or_") != std::string::npos)   side = "OR(#eta_{1}, #eta_{2})";
        else if (metric_name.find("and_") != std::string::npos)  side = "AND(#eta_{1}, #eta_{2})";

        if (metric_name.find("external") != std::string::npos)   trigger = "External Trigger";
        else if (metric_name.find("rpc") != std::string::npos)   trigger = "RPC Coincidence";

        // Build a subtitle string
        // This applies only for those metrics that depend on the track reconstruction step
        std::smatch match;
        if (std::regex_search(metric_name, match, std::regex("^(track_)?(avg_tot|avg_multiplicity|eff)_"))) {
            out_title += match[1].matched ? "After Track Reco: " : "Before Track Reco: ";
        }

        if (std::regex_search(metric_name, match, std::regex("^(track_)?(avg_tot|avg_multiplicity|noise_rate_strips)_"))) {
            out_title += "Layer " + matchLabels(metric_name, "layer(\\d+)") + ": ";
        }

        if (!side.empty()) out_title += side;
        if (!trigger.empty()) {
            if (!out_title.empty()) out_title += ": ";
            out_title += trigger;
        }

        return std::make_tuple(out_title, out_xaxis, out_yaxis, legend_entries);
    }

    enum class AxisType { X, Y, Z };

    void setRange(TObject* obj, TAxis* axis, AxisType axis_type,
        double default_min, double default_max, double x_floor_limit = -1e9) {
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
                    if (x_val <= x_floor_limit) continue;

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
                    } else { // Y axis
                        val_low = val_high = gr->GetY()[i];
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
                    if (x_center <= x_floor_limit) continue;

                    // For 1D histograms, Y and Z loops run exactly once (GetNbinsY/Z return 1)
                    for (int y = 1; y <= n_bins_y; ++y) {
                        for (int z = 1; z <= n_bins_z; ++z) {
                            int global_bin = h->GetBin(x, y, z);
                            double content = h->GetBinContent(global_bin);
                            double error = h->GetBinError(global_bin);

                            // Skip perfectly empty bins so they don't drag the Y-axis to 0
                            if (content == 0 && error == 0) continue; 

                            double val_low = 0.0, val_high = 0.0;
                            if (axis_type == AxisType::X) {
                                val_low = h->GetXaxis()->GetBinLowEdge(x);
                                val_high = h->GetXaxis()->GetBinUpEdge(x);
                            } else {
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
            double safety_buffer = (true_max > true_min) ? (true_max - true_min) * 0.05 : 50.0;

            double dynamic_min = true_min - safety_buffer;
            double dynamic_max = true_max + safety_buffer;

            // For Y-axis, clamp to the absolute floor (e.g., 0.0) if it shouldn't dip negative
            if (axis_type == AxisType::Y && dynamic_min < default_min) {
                dynamic_min = default_min; 
            }

            axis->SetLimits(dynamic_min, dynamic_max);
            axis->SetRangeUser(dynamic_min, dynamic_max);
        } else {
            // Fallback to absolute defaults if no valid points pass the threshold
            axis->SetLimits(default_min, default_max);
            axis->SetRangeUser(default_min, default_max);
        }
    }

    // --------------------------------------------------------------------------------------
    namespace ATLASStyler {

        const std::vector<Color_t> ATLAS_PALETTE = {
            static_cast<Color_t>(TColor::GetColor("#144d92")),    // Sharp Blue
            static_cast<Color_t>(TColor::GetColor("#CF4446")),    // Deep Red
            static_cast<Color_t>(TColor::GetColor("#1a8f3f")),    // Dark Green
            static_cast<Color_t>(TColor::GetColor("#e28843"))     // Burnt Orange
        };

        std::vector<TLatex*> drawATLASLabel(float ndc_x, float ndc_y,
            const std::string& status, short alignment = 11) {
            std::vector<TLatex*> drawn_objects;
            if (!gPad) return drawn_objects;

            TLatex* l_atlas = new TLatex();
            l_atlas->SetNDC();
            l_atlas->SetTextColor(kBlack);
            l_atlas->SetTextSize(0.04);
            l_atlas->SetTextAlign(alignment);

            if (status.empty()) {
                l_atlas->SetTextFont(72);
                l_atlas->SetTextAlign(alignment);
                l_atlas->DrawLatex(ndc_x, ndc_y, "ATLAS");
                drawn_objects.push_back(l_atlas);
            } else {
                TLatex* l_status = new TLatex();
                l_status->SetNDC();
                l_status->SetTextColor(kBlack);
                l_status->SetTextSize(0.04);
                l_status->SetTextAlign(alignment);

                if (alignment >= 30 && alignment <= 33) { // Right-aligned family
                    l_status->SetTextFont(42);
                    TLatex* s_drawn = l_status->DrawLatex(ndc_x, ndc_y, status.c_str());

                    double status_width = 0.12;
                    if (s_drawn) {
                        UInt_t w = 0, h = 0;
                        s_drawn->GetBoundingBox(w, h);
                        status_width = static_cast<double>(w) / gPad->GetWw();
                    }

                    l_atlas->SetTextFont(72);
                    TLatex* a_drawn = l_atlas->DrawLatex(ndc_x - status_width - 0.01, ndc_y, "ATLAS");

                    drawn_objects.push_back(a_drawn ? a_drawn : l_atlas);
                    drawn_objects.push_back(s_drawn ? s_drawn : l_status);
                } else { // Left-aligned family
                    l_atlas->SetTextFont(72);
                    TLatex* a_drawn = l_atlas->DrawLatex(ndc_x, ndc_y, "ATLAS");

                    double atlas_width = 0.12;
                    if (a_drawn) {
                        UInt_t w = 0, h = 0;
                        a_drawn->GetBoundingBox(w, h);
                        atlas_width = static_cast<double>(w) / gPad->GetWw();
                    }

                    l_status->SetTextFont(42);
                    TLatex* s_drawn = l_status->DrawLatex(ndc_x + atlas_width + 0.01, ndc_y, status.c_str());

                    drawn_objects.push_back(a_drawn ? a_drawn : l_atlas);
                    drawn_objects.push_back(s_drawn ? s_drawn : l_status);
                }
            }

            return drawn_objects;
        }

        TLatex* drawPlotTitle(TObject* obj, float ndc_x, float ndc_y, short alignment = 11) {
            if (!obj || !gPad) return nullptr;

            std::string titleStr = obj->GetTitle();
            if (titleStr.empty()) return nullptr;

            TLatex* t = new TLatex();
            t->SetNDC();
            t->SetTextFont(42);
            t->SetTextSize(0.04);
            t->SetTextColor(kBlack);
            t->SetTextAlign(alignment);

            TLatex* title_drawn = t->DrawLatex(ndc_x, ndc_y, titleStr.c_str());

            return title_drawn ? title_drawn : t; 
        }

        TPaveText* drawATLASHeaderBlock(
            double ndc_x, double ndc_y,
            const std::string& status = "",
            const std::string& title = "",
            short alignment = 33,
            Color_t fillColor = kWhite, double fillAlpha = 0.70,
            Color_t borderColor = kBlack, int borderWidth = 1,
            double innerPadding = 0.01)
        {
            if (!gPad) return nullptr;

            // Construct text strings
            std::string line1 = "#bf{#it{ATLAS}}";
            if (!status.empty()) {
                line1 += " " + status;
            }

            // Measure actual text width dynamically on gPad using a dummy TLatex
            auto getWidthNDC = [&](const std::string& txt) -> double {
                if (txt.empty()) return 0.0;

                TLatex measure;
                measure.SetNDC();
                measure.SetTextFont(42);
                measure.SetTextSize(0.035);
                measure.SetText(0, 0, txt.c_str());

                UInt_t w = 0, h = 0;
                measure.GetBoundingBox(w, h);
                return static_cast<double>(w) / gPad->GetWw();
            };

            double w1 = getWidthNDC(line1);
            double w2 = getWidthNDC(title);
            double max_text_width = std::max(w1, w2);

            int num_lines = title.empty() ? 1 : 2;
            double line_height = 0.038;

            // Calculate total box dimensions including padding
            double box_width  = max_text_width + (2.0 * innerPadding);
            double box_height = (num_lines * line_height) + (2.0 * innerPadding);

            // Compute (x1, y1, x2, y2) relative to (ndc_x, ndc_y) using the alignment code
            double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
            int h_align = alignment / 10; // 1 = Left, 2 = Center, 3 = Right
            int v_align = alignment % 10; // 1 = Bottom, 2 = Middle, 3 = Top

            if (h_align == 3)      { x2 = ndc_x; x1 = ndc_x - box_width; }
            else if (h_align == 2) { x1 = ndc_x - (box_width / 2.0); x2 = ndc_x + (box_width / 2.0); }
            else                   { x1 = ndc_x; x2 = ndc_x + box_width; }

            if (v_align == 3)      { y2 = ndc_y; y1 = ndc_y - box_height; }
            else if (v_align == 2) { y1 = ndc_y - (box_height / 2.0); y2 = ndc_y + (box_height / 2.0); }
            else                   { y1 = ndc_y; y2 = ndc_y + box_height; }

            // Instantiate TPaveText with auto-fitted coordinates
            TPaveText* pave = new TPaveText(x1, y1, x2, y2, "NDC");

            pave->SetMargin(0.02);
            pave->SetBorderSize(borderWidth);
            pave->SetLineColor(borderColor);
            pave->SetLineWidth(borderWidth);
            pave->SetLineStyle(1);
            pave->SetShadowColor(0);

            if (fillAlpha < 1.0) {
                pave->SetFillColorAlpha(fillColor, fillAlpha);
            } else {
                pave->SetFillColor(fillColor);
            }

            pave->SetTextAlign(alignment);
            pave->SetTextFont(42);
            pave->SetTextSize(0.035);

            // Add lines
            TText* t1 = pave->AddText(line1.c_str());
            if (t1) t1->SetTextColor(kBlack);

            if (!title.empty()) {
                TText* t2 = pave->AddText(title.c_str());
                if (t2) t2->SetTextColor(kBlack);
            }

            pave->Draw();
            return pave;
        }

        TLegend* drawATLASLegend(TObject* obj, const std::vector<std::string>& legend_entries,
            float ndc_x, float ndc_y, short alignment) {

            if (!obj || !gPad) return nullptr;

            TList* items_list = nullptr;
            bool is_stack = false;

            // Check container type
            if (auto* mg = dynamic_cast<TMultiGraph*>(obj)) {
                items_list = mg->GetListOfGraphs();
            } else if (auto* stack = dynamic_cast<THStack*>(obj)) {
                items_list = stack->GetHists();
                is_stack = true;
            }

            if (!items_list || items_list->GetSize() == 0) return nullptr;

            // Dimensions of the legend box
            double leg_width = 0.25; 
            double entry_height = 0.04; 
            double leg_height = items_list->GetSize() * entry_height;

            double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;

            // Parse the 2-digit alignment code
            int h_align = alignment / 10;
            int v_align = alignment % 10;

            // Horizontal anchoring
            if (h_align == 3) {
                x2 = ndc_x; x1 = ndc_x - leg_width;
            } else if (h_align == 2) {
                x1 = ndc_x - (leg_width / 2.0); x2 = ndc_x + (leg_width / 2.0);
            } else {
                x1 = ndc_x; x2 = ndc_x + leg_width;
            }

            // Vertical anchoring
            if (v_align == 3) {
                y2 = ndc_y; y1 = ndc_y - leg_height;
            } else if (v_align == 2) {
                y1 = ndc_y - (leg_height / 2.0); y2 = ndc_y + (leg_height / 2.0);
            } else {
                y1 = ndc_y; y2 = ndc_y + leg_height;
            }

            TLegend* leg = new TLegend(x1, y1, x2, y2);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.04);

            TIter next_item(items_list);
            TObject* child = nullptr;
            size_t idx = 0;

            while ((child = next_item())) {
                // Fallback to title if the vector is empty or out of bounds
                std::string label = child->GetTitle();

                // Inject the dynamically compiled label if available
                if (idx < legend_entries.size() && !legend_entries[idx].empty()) {
                    label = legend_entries[idx];
                }

                std::string draw_option = is_stack ? "f" : "pe";
                leg->AddEntry(child, label.c_str(), draw_option.c_str());

                idx++;
            }

            leg->Draw();
            return leg;
        }
        // ----------------------------------------------------------------------------------
        namespace {

            struct SpacingSettings {
                double margin;      // Suggested pad margin (left for Y, right for Z, bottom for X)
                double titleOffset; // Suggested axis title offset
            };

            int getLabelWidthDigits(double val) {
                double absVal = std::abs(val);
                if (absVal == 0.0) return 1;

                if (absVal >= 1.0) {
                    return static_cast<int>(std::floor(std::log10(absVal))) + 1;
                } else {
                    int decimalPlaces = static_cast<int>(-std::floor(std::log10(absVal)));
                    return decimalPlaces + 2;
                }
            }

            // Computes digits and returns spacing settings
            SpacingSettings calculateAxisSpacing(double maxVal, AxisType type) {
                int digits = getLabelWidthDigits(maxVal);

                SpacingSettings settings;

                // Map digits to margins and offsets based on axis orientation
                if (type == AxisType::Z) { // Colorbar (Right margin & Z-title offset)
                    if (digits > 5)       { settings.margin = 0.18; settings.titleOffset = 1.45; }
                    else if (digits >= 4) { settings.margin = 0.17; settings.titleOffset = 1.15; }
                    else if (digits == 3) { settings.margin = 0.16; settings.titleOffset = 1.00; }
                    else if (digits == 2) { settings.margin = 0.15; settings.titleOffset = 0.85; }
                    else                  { settings.margin = 0.14; settings.titleOffset = 0.70; }
                } else if (type == AxisType::Y) { // Y-Axis (Left margin & Y-title offset)
                    if (digits > 5)       { settings.margin = 0.18; settings.titleOffset = 1.40; }
                    else if (digits >= 4) { settings.margin = 0.16; settings.titleOffset = 1.25; }
                    else if (digits == 3) { settings.margin = 0.14; settings.titleOffset = 1.10; }
                    else if (digits == 2) { settings.margin = 0.12; settings.titleOffset = 0.95; }
                    else                  { settings.margin = 0.10; settings.titleOffset = 0.80; }
                } else { // X-Axis (Bottom margin & X-title offset)
                    settings.margin = 0.14;
                    settings.titleOffset = 1.10;
                }

                return settings;
            }

            void setMargins(TObject* obj, TPad* pad, AxisType axisType) {
                if (!obj || !pad) return;

                // Set default margins to be small
                pad->SetTopMargin(0.05);
                pad->SetRightMargin(0.05);

                double maxVal = 0.0;
                TAxis* axis = nullptr;

                // Extract axis and maximum value depending on object type
                if (auto h2 = dynamic_cast<TH2*>(obj)) {
                    if (axisType == AxisType::X)      { axis = h2->GetXaxis(); maxVal = h2->GetXaxis()->GetXmax(); }
                    else if (axisType == AxisType::Y) { axis = h2->GetYaxis(); maxVal = h2->GetYaxis()->GetXmax(); }
                    else if (axisType == AxisType::Z) { axis = h2->GetZaxis(); maxVal = h2->GetMaximum(); }
                } else if (auto h1 = dynamic_cast<TH1*>(obj)) {
                    if (axisType == AxisType::X)      { axis = h1->GetXaxis(); maxVal = h1->GetXaxis()->GetXmax(); }
                    else if (axisType == AxisType::Y) { axis = h1->GetYaxis(); maxVal = h1->GetMaximum(); }
                } else if (auto stack = dynamic_cast<THStack*>(obj)) {
                    if (axisType == AxisType::X)      { axis = stack->GetXaxis(); if (axis) maxVal = axis->GetXmax(); }
                    else if (axisType == AxisType::Y) { axis = stack->GetYaxis(); maxVal = stack->GetMaximum(); }
                } else if (auto eg = dynamic_cast<TGraphAsymmErrors*>(obj)) {
                    if (axisType == AxisType::X)      { axis = eg->GetXaxis(); maxVal = eg->GetXaxis()->GetXmax(); }
                    else if (axisType == AxisType::Y) { axis = eg->GetYaxis(); maxVal = eg->GetYaxis()->GetXmax(); }
                } else if (auto g = dynamic_cast<TGraph*>(obj)) {
                    if (axisType == AxisType::X)      { axis = g->GetXaxis(); maxVal = g->GetXaxis()->GetXmax(); }
                    else if (axisType == AxisType::Y) { axis = g->GetYaxis(); maxVal = g->GetYaxis()->GetXmax(); }
                } else if (auto mg = dynamic_cast<TMultiGraph*>(obj)) {
                    if (axisType == AxisType::X)      { axis = mg->GetXaxis(); maxVal = mg->GetXaxis()->GetXmax(); }
                    else if (axisType == AxisType::Y) { axis = mg->GetYaxis(); maxVal = mg->GetYaxis()->GetXmax(); }
                }

                if (!axis) return;

                // Calculate settings
                SpacingSettings spacing = calculateAxisSpacing(maxVal, axisType);

                // Apply pad margins depending on axis orientation
                if (axisType == AxisType::Z)      pad->SetRightMargin(spacing.margin);
                else if (axisType == AxisType::Y) pad->SetLeftMargin(spacing.margin);
                else if (axisType == AxisType::X) pad->SetBottomMargin(spacing.margin);

                // Set axis title offset
                axis->SetTitleOffset(spacing.titleOffset);

                pad->Modified();
                pad->Update();
            }

            // Adjusts the color bar of a 2D histogram dynamically based on the maximum value
            void adjustDynamicCB(TH2* h2, TPad* pad) {
                if (!h2 || !pad) return;

                // Delegate dynamic right margin & Z-title offset calculation to estimateMargins
                setMargins(h2, pad, AxisType::Z);

                // Find the colorbar object
                TObject* paletteObj = h2->GetListOfFunctions() ? h2->GetListOfFunctions()->FindObject("palette") : nullptr;
                if (paletteObj) {
                    // Find where the plot frame ends
                    double rightMargin = pad->GetRightMargin();
                    double plotFrameRightX = 1.0 - rightMargin;

                    // Anchor the colorbar a constant distance away from the plot frame
                    double gap = 0.015;
                    double barWidth = 0.035;

                    double x1 = plotFrameRightX + gap;
                    double x2 = plotFrameRightX + gap + barWidth;
                    double y1 = pad->GetBottomMargin();
                    double y2 = 1.0 - pad->GetTopMargin();

                    // Apply coordinates dynamically
                    char paramStr[64];

                    snprintf(paramStr, sizeof(paramStr), "%f", x1);
                    paletteObj->Execute("SetX1NDC", paramStr);

                    snprintf(paramStr, sizeof(paramStr), "%f", x2);
                    paletteObj->Execute("SetX2NDC", paramStr);

                    snprintf(paramStr, sizeof(paramStr), "%f", y1);
                    paletteObj->Execute("SetY1NDC", paramStr);

                    snprintf(paramStr, sizeof(paramStr), "%f", y2);
                    paletteObj->Execute("SetY2NDC", paramStr);
                }

                pad->Modified();
                pad->Update();
            }

            std::array<TAxis*, 3> extractAxis(TObject* obj, TPad* pad) {
                std::array<TAxis*, 3> axes = {nullptr, nullptr, nullptr};

                // Safely extract axes depending on the object type
                if (auto h2 = dynamic_cast<TH2*>(obj)) {
                    axes[0] = h2->GetXaxis();
                    axes[1] = h2->GetYaxis();
                    axes[2] = h2->GetZaxis();
                } else if (auto h = dynamic_cast<TH1*>(obj)) {
                    axes[0] = h->GetXaxis();
                    axes[1] = h->GetYaxis();
                    h->SetMarkerStyle(20);
                    h->SetMarkerSize(1.0);
                    h->SetLineColor(kBlack);
                    h->SetLineWidth(2);
                } else if (auto stack = dynamic_cast<THStack*>(obj)) {
                    TH1* frame = stack->GetHistogram();
                    if (frame) {
                        axes[0] = frame->GetXaxis();
                        axes[1] = frame->GetYaxis();
                    }
                } else if (auto mg = dynamic_cast<TMultiGraph*>(obj)) {
                    axes[0] = mg->GetXaxis();
                    axes[1] = mg->GetYaxis();
                    TList* graph_list = mg->GetListOfGraphs();
                    if (!graph_list) return axes;
                }
                else if (auto g = dynamic_cast<TGraph*>(obj)) {
                    axes[0] = g->GetXaxis();
                    axes[1] = g->GetYaxis();
                    g->SetMarkerStyle(20);
                    g->SetMarkerSize(1.0);
                    g->SetLineColor(kBlack);
                    g->SetLineWidth(2);
                }

                pad->Modified();
                pad->Update();

                return axes;
            }

            void styleAxis(TAxis* axis, TPad* pad, TObject* parentObj, AxisType type,
                float labelSize = 0.04f, float titleSize = 0.05f, int ndivisions = 510) {
                if (!axis || !pad) return;

                axis->SetLabelFont(42);
                axis->SetTitleFont(42);
                axis->SetLabelSize(labelSize);
                axis->SetTitleSize(titleSize);
                axis->SetNdivisions(ndivisions, kTRUE);

                // Dynamically set offset & pad margins if parent object is provided
                if (parentObj) setMargins(parentObj, pad, type);

                pad->Modified();
                pad->Update();
            }
        }   // anonymous namespace
        // ----------------------------------------------------------------------------------

        void applyATLASStyle(TObject* obj, TPad* pad) {

            gStyle->SetOptTitle(0);
            gStyle->SetOptStat(0);

            // Custom dash lines: dash/gap/dot/gap
            gStyle->SetLineStyleString(11, "48 24");

            if (gPad) {
                gPad->SetTickx(1); // 1 = Draw ticks on top side
                gPad->SetTicky(1); // 1 = Draw ticks on right side
            }

            // Extract axes for styling
            auto [xAxis, yAxis, zAxis] = extractAxis(obj, pad);

            // Apply ATLAS standard font rules to all axes
            if (xAxis) styleAxis(xAxis, pad, obj, AxisType::X, 0.04, 0.05, 510);
            if (yAxis) styleAxis(yAxis, pad, obj, AxisType::Y, 0.04, 0.05, 510);

            // For 2D histograms, apply dynamic right margin and Z-axis title offset
            if (zAxis) styleAxis(zAxis, pad, nullptr, AxisType::Z, 0.04, 0.05, 510);

            // Adjust the colorbar for 2D histograms
            if (auto h2 = dynamic_cast<TH2*>(obj)) adjustDynamicCB(h2, pad);

            // Set colormap and palette for 2D histograms
            if (zAxis) gStyle->SetPalette(kBird);

            pad->Modified();
            pad->Update();
        }

    }   // namespace ATLASStyler
    // --------------------------------------------------------------------------------------

    namespace Objects {

        void line(TVirtualPad* canvas, float x_start, float x_end, float y_start, float y_end,
            Color_t color = kBlack, int line_style = 1, int line_width = 1) {
            if (!canvas) return;

            canvas->cd();
            TLine* line = new TLine(x_start, y_start, x_end, y_end);
            line->SetLineColor(color);
            line->SetLineStyle(line_style);
            line->SetLineWidth(line_width);
            line->Draw("same");
        }

        void box(TVirtualPad* canvas, float x_start, float x_end, float y_start, float y_end,
                Color_t color = kWhite, float alpha = 1.0,
                int line_width = 1, int line_style = 1, Color_t line_color = kBlack, float line_alpha = 1.0)
        {
            if (!canvas) return;
            canvas->cd();

            if (x_start >= x_end || y_start >= y_end) {
                std::cerr << "Error: Invalid coordinates for box. Ensure x_start < x_end and y_start < y_end." << std::endl;
                return;
            }

            TPave* box = new TPave(x_start, y_start, x_end, y_end, line_width, "NDC");
            box->SetLineWidth(line_width);
            box->SetLineStyle(line_style);
            if (alpha < 1.0f) {
                box->SetFillColorAlpha(color, alpha);
                box->SetLineColorAlpha(line_color, line_alpha);
            } else {
                box->SetFillColor(color);
                box->SetLineColor(line_color);
            }

            box->Draw("same");

            canvas->Modified();
            canvas->Update();
        }

        void hatchedRegion(TVirtualPad* canvas,
            float x_start, float x_end, float y_start, float y_end,
            int pattern = 3245, int line_width = 1,
            Color_t color = kAzure + 2, float alpha = 0.5)
        {
            canvas->cd();

            // Check if the coordinates are valid
            if (x_start >= x_end || y_start >= y_end) {
                std::cerr << "Error: Invalid coordinates for hatched region. Ensure x_start < x_end and y_start < y_end." << std::endl;
                return;
            }

            TBox* box = new TBox(x_start, y_start, x_end, y_end);
            box->SetFillStyle(pattern);

            if (alpha < 1.0) {
                box->SetFillColorAlpha(color, alpha);
                box->SetLineColorAlpha(color, alpha);
            } else {
                box->SetFillColor(color);
                box->SetLineColor(color);
            }

            box->SetLineWidth(line_width);

            box->Draw("same");

            canvas->Modified();
            canvas->Update();
        }

        std::pair<float, float> getTextSizeNDC(TLatex* latex_obj) {
            if (!latex_obj || !gPad) return {0.0, 0.0};
            UInt_t w = 0, h = 0;
            latex_obj->GetBoundingBox(w, h);
            return {static_cast<double>(w) / gPad->GetWw(), static_cast<double>(h) / gPad->GetWh()};
        }

    }   // namespace Objects
    // --------------------------------------------------------------------------------------

    using namespace ATLASStyler;

    void styleEfficiencyVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        // Extract title and axis labels from the object's title string
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        obj->Draw("AP0Z");

        // Set axis ranges and labels
        if (mg && mg->GetHistogram()) {

            // Setup X Axis
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, 0.0, 10e3, 4500.0);
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

    void styleAvgToTVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        // Differentiate between a 24-strip plot and a 3-layer plot
        bool is_strip_plot = (mg && mg->GetListOfGraphs() && mg->GetListOfGraphs()->GetSize() > 10);

        if (is_strip_plot) gStyle->SetPalette(kViridis);

        obj->Draw("APZ");

        if (mg && mg->GetHistogram()) {
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, 0.0, 10e3, 4500.0);
                xAxis->SetTitle(x_label.c_str());
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, 0.0, 1.0);
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
                setRange(mg, xAxis, AxisType::X, 0.0, 10e3, 4500.0);
                xAxis->SetTitle(x_label.c_str());
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, 0.8, 2.0);
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

        double ndc_x0 = pad1->GetLeftMargin();
        double ndc_y0 = 1.0 - pad1->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.92,
            ndc_y0 - 0.05,            // Coordinates for the header box
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

        auto enforceIntegerMinorTicks = [](TAxis* axis) {
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

    static const std::vector<std::pair<PlotCategory, StylerFnPtr>> styler_map = {
        {PlotCategory::EfficiencyVsHV,          &styleEfficiencyVsHV},
        {PlotCategory::MeanClusterSizeVsHV,     &styleEfficiencyVsHV},
        {PlotCategory::NoiseRateVsHV,           &styleEfficiencyVsHV},
        {PlotCategory::AvgToTVsHV,              &styleAvgToTVsHV},
        {PlotCategory::AvgMultVsHV,             &styleAvgMulVsHV},
        {PlotCategory::StripDistribution,       &styleStripDistribution},
        {PlotCategory::ToTDistribution,         &styleToTDistribution},
        {PlotCategory::ToTCombinedDistribution, &styleToTCombinedDistribution},
        {PlotCategory::Default,                 &styleDefaultPlot}
    };

    StylerFnPtr getCustomStyler(PlotCategory category) {
        for (const auto& [cat, fn_ptr] : styler_map) {
            if (cat == category) return fn_ptr;
        }
        return nullptr;
    }
}   // namespace PlotStyler
}   // namespace PlotterHelpers