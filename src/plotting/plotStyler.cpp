#include <cmath>
#include <utility>

#include <iostream>

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
        {"h1d_strip_eta",    TH1::Class(),                  PlotCategory::StripDistribution},
        {"h1d_tot_eta",      TH1::Class(),                  PlotCategory::ToTDistribution},
        {"h1d_tot",          THStack::Class(),              PlotCategory::ToTCombinedDistribution},
        {"track_eff",        TGraphAsymmErrors::Class(),    PlotCategory::Efficiency},
        {"eff",              TGraphAsymmErrors::Class(),    PlotCategory::Efficiency},
        {"track_eff",        TMultiGraph::Class(),          PlotCategory::EfficiencyVsHV},
        {"eff",              TMultiGraph::Class(),          PlotCategory::EfficiencyVsHV},
        {"avg_cluster_size", TMultiGraph::Class(),          PlotCategory::MeanClusterSizeVsHV},
        {"noise_rate",       TMultiGraph::Class(),          PlotCategory::NoiseRateVsHV}
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

    std::tuple<std::string, std::string, std::string> compilePlotLabels(const std::string& metric_name) {
        std::string out_title = "", out_xaxis = "", out_yaxis = "";
        std::string side = "", trigger = "";

        if (metric_name.rfind("track_eff_", 0) == 0 || metric_name.rfind("eff_", 0) == 0) {
            out_xaxis = "High Voltage [V]"; out_yaxis = "Efficiency";
        } else if (metric_name.rfind("avg_cluster_size_", 0) == 0) {
            out_xaxis = "High Voltage [V]"; out_yaxis = "#LTCluster Size#GT [Hits]";
        } else if (metric_name.rfind("noise_rate_", 0) == 0) {
            out_xaxis = "High Voltage [V]"; out_yaxis = "Noise Rate [Hz/cm^{2}]";
        }

        if (metric_name.find("eta1") != std::string::npos)       side = "#eta_{1} Side";
        else if (metric_name.find("eta2") != std::string::npos)  side = "#eta_{2} Side";
        else if (metric_name.find("or_") != std::string::npos)   side = "OR(#eta_{1}, #eta_{2})";
        else if (metric_name.find("and_") != std::string::npos)  side = "AND(#eta_{1}, #eta_{2})";

        if (metric_name.find("external") != std::string::npos)   trigger = "External Trigger";
        else if (metric_name.find("rpc") != std::string::npos)   trigger = "RPC Coincidence";

        if (!side.empty()) out_title += side;
        if (!trigger.empty()) {
            if (!out_title.empty()) out_title += ": ";
            out_title += trigger;
        }
        return std::make_tuple(out_title, out_xaxis, out_yaxis);
    }

    // --------------------------------------------------------------------------------------
    namespace ATLASStyler {

        const std::vector<Color_t> ATLAS_PALETTE = {
            static_cast<Color_t>(TColor::GetColor("#144d92")),    // Sharp Blue
            static_cast<Color_t>(TColor::GetColor("#CF4446")),    // Deep Red
            static_cast<Color_t>(TColor::GetColor("#1a8f3f")),    // Dark Green
            static_cast<Color_t>(TColor::GetColor("#e28843"))     // Burnt Orange
        };

        std::vector<TLatex*> drawATLASLabel(float ndc_x, float ndc_y, const std::string& status, short alignment = 11) {
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

        TPaveText* drawATLASHeaderBlock(double ndc_x, double ndc_y,
                                        const std::string& status = "",
                                        const std::string& title = "",
                                        short alignment = 33,
                                        Color_t fillColor = kWhite, float fillAlpha = 0.70f,
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

            if (fillAlpha < 1.0f) {
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

        TLegend* drawATLASLegend(TObject* obj, float ndc_x, float ndc_y, short alignment = 11) {
            if (!obj || !gPad) return nullptr;

            TList* items_list = nullptr;
            bool is_graph = false;
            bool is_stack = false;

            // Check container type
            if (auto* mg = dynamic_cast<TMultiGraph*>(obj)) {
                items_list = mg->GetListOfGraphs();
                is_graph = true;
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

            // Parse the 2-digit alignment code:
            // First digit:  1 = Left, 2 = Center, 3 = Right
            // Second digit: 1 = Bottom, 2 = Middle, 3 = Top
            int h_align = alignment / 10;
            int v_align = alignment % 10;

            // Horizontal anchoring
            if (h_align == 3) {        // Right-aligned: ndc_x is the RIGHT edge
                x2 = ndc_x;
                x1 = ndc_x - leg_width;
            } else if (h_align == 2) { // Centered: ndc_x is the CENTER
                x1 = ndc_x - (leg_width / 2.0);
                x2 = ndc_x + (leg_width / 2.0);
            } else {                   // Left-aligned (1 or default): ndc_x is the LEFT edge
                x1 = ndc_x;
                x2 = ndc_x + leg_width;
            }

            // Vertical anchoring
            if (v_align == 3) {        // Top-aligned: ndc_y is the TOP edge
                y2 = ndc_y;
                y1 = ndc_y - leg_height;
            } else if (v_align == 2) { // Middle-aligned: ndc_y is the CENTER
                y1 = ndc_y - (leg_height / 2.0);
                y2 = ndc_y + (leg_height / 2.0);
            } else {                   // Bottom-aligned (1 or default): ndc_y is the BOTTOM edge
                y1 = ndc_y;
                y2 = ndc_y + leg_height;
            }

            TLegend* leg = new TLegend(x1, y1, x2, y2);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.04);

            TIter next_item(items_list);
            TObject* child = nullptr;

            while ((child = next_item())) {
                std::string name = child->GetName();
                std::string label = child->GetTitle(); // Use configured title by default

                if (is_graph) {
                    if (name.find("layer0") != std::string::npos) label = "Layer 0";
                    else if (name.find("layer1") != std::string::npos) label = "Layer 1";
                    else if (name.find("layer2") != std::string::npos) label = "Layer 2";

                    leg->AddEntry(child, label.c_str(), "pe");
                } 
                else if (is_stack) {
                    if (name.find("tot_eta1") != std::string::npos) label = "#eta1 Side";
                    else if (name.find("tot_eta2") != std::string::npos) label = "#eta2 Side";
                    
                    leg->AddEntry(child, label.c_str(), "f");
                }
            }

            leg->Draw();
            return leg;
        }

        void adjustDynamicCB(TH2* h2, TPad* pad) {
            if (!h2 || !pad) return;

            double maxVal = h2->GetMaximum();
            int digits = 1;
            if (maxVal > 0) {
                digits = static_cast<int>(std::log10(maxVal)) + 1;
            }

            double rightMargin = 0.15;
            double titleOffset = 0.95;

            if (digits > 5) {
                rightMargin = 0.19;
                titleOffset = 1.35;
            } else if (digits >= 4) {
                rightMargin = 0.17; 
                titleOffset = 1.10;
            }

            // Apply directly to the pad pointer before any Draw() execution happens
            pad->SetRightMargin(rightMargin);
            
            TAxis* zAxis = h2->GetZaxis();
            if (zAxis) {
                zAxis->SetTitleOffset(titleOffset);
            }
        }

        void applyATLASStyle(TObject* obj, TPad* pad) {
            TAxis* xAxis = nullptr;
            TAxis* yAxis = nullptr;
            TAxis* zAxis = nullptr;

            gStyle->SetOptTitle(0);
            gStyle->SetOptStat(0);

            // Custom dash lines: dash/gap/dot/gap
            gStyle->SetLineStyleString(11, "48 24");

            if (gPad) {
                gPad->SetTickx(1); // 1 = Draw ticks on top side
                gPad->SetTicky(1); // 1 = Draw ticks on right side
            }

            // Safely extract axes depending on the object type
            if (auto h2 = dynamic_cast<TH2*>(obj)) {
                xAxis = h2->GetXaxis();
                yAxis = h2->GetYaxis();
                zAxis = h2->GetZaxis();
                adjustDynamicCB(h2, pad); // Adjust color bar dynamically based on max value
            } else if (auto h = dynamic_cast<TH1*>(obj)) {
                xAxis = h->GetXaxis();
                yAxis = h->GetYaxis();
                h->SetMarkerStyle(20);
                h->SetMarkerSize(1.0);
                h->SetLineColor(kBlack);
                h->SetLineWidth(2);
            } else if (auto stack = dynamic_cast<THStack*>(obj)) {
                TH1* frame = stack->GetHistogram();
                if (frame) {
                    xAxis = frame->GetXaxis();
                    yAxis = frame->GetYaxis();
                }
            } else if (auto mg = dynamic_cast<TMultiGraph*>(obj)) {
                xAxis = mg->GetXaxis();
                yAxis = mg->GetYaxis();
                TList* graph_list = mg->GetListOfGraphs();
                if (!graph_list) return;

                TIter next_graph(graph_list);
                TGraph* sub_graph = nullptr;
                int color_index = 0;

                while ((sub_graph = static_cast<TGraph*>(next_graph()))) {
                    sub_graph->SetMarkerStyle(5);
                    sub_graph->SetMarkerSize(1.0);

                    std::string g_name = sub_graph->GetName();
                    int layer_id = 0; // Default fallback
                    
                    if (g_name.find("layer0") != std::string::npos) layer_id = 0;
                    else if (g_name.find("layer1") != std::string::npos) layer_id = 1;
                    else if (g_name.find("layer2") != std::string::npos) layer_id = 2;

                    Color_t assigned_color = ATLAS_PALETTE[0]; // Default to kBlack
                    if (layer_id >= 0 && layer_id < static_cast<int>(ATLAS_PALETTE.size())) {
                        assigned_color = ATLAS_PALETTE[layer_id];
                    }

                    sub_graph->SetMarkerColor(assigned_color);
                    sub_graph->SetLineColor(assigned_color);

                    color_index++;
                }
            }
            else if (auto g = dynamic_cast<TGraph*>(obj)) {
                xAxis = g->GetXaxis();
                yAxis = g->GetYaxis();
                g->SetMarkerStyle(20);
                g->SetMarkerSize(1.0);
                g->SetLineColor(kBlack);
                g->SetLineWidth(2);
            }

            // Set colormap and palette for 2D histograms
            gStyle->SetPalette(kBird);

            // Apply ATLAS standard font rules (Font 42 = Helvetica, Font 43 = Precise Pixel size) to all axes
            if (xAxis && yAxis) {
                xAxis->SetLabelFont(42);
                xAxis->SetTitleFont(42);
                xAxis->SetLabelSize(0.04);
                xAxis->SetTitleSize(0.05);
                xAxis->SetTitleOffset(1.1);
                xAxis->SetNdivisions(510, kTRUE);

                yAxis->SetLabelFont(42);
                yAxis->SetTitleFont(42);
                yAxis->SetLabelSize(0.04);
                yAxis->SetTitleSize(0.05);
                yAxis->SetTitleOffset(1.3);
                yAxis->SetNdivisions(510, kTRUE);
            }

            if (zAxis) {
                zAxis->SetLabelFont(42);
                zAxis->SetTitleFont(42);
                zAxis->SetLabelSize(0.04);
                zAxis->SetTitleSize(0.05);
            }
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
        // Margins
        canvas->SetLeftMargin(0.15);
        canvas->SetRightMargin(0.05);
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.14);

        // Extract title and axis labels from the object's title string
        auto [title, x_label, y_label] = compilePlotLabels(obj->GetTitle());

        // Draw first to generate internal axis frame
        obj->Draw("APZ");

        auto mg = dynamic_cast<TMultiGraph*>(obj);

        // Set axis ranges and labels
        if (mg && mg->GetHistogram()) {
            double floor_limit = 4500.0; // Keep floor_limit in parent scope for both axes

            TAxis* xAxis = mg->GetHistogram()->GetXaxis();
            if (xAxis) {
                double dynamic_lower_bound = floor_limit;
                double dynamic_upper_bound = floor_limit;

                if (mg->GetListOfGraphs()) {
                    double true_min_x = INT_MAX;
                    double true_max_x = INT_MIN;

                    TIter next(mg->GetListOfGraphs());
                    TObject* gr_obj = nullptr;

                    while ((gr_obj = next())) {
                        auto* gr = dynamic_cast<TGraph*>(gr_obj);
                        if (gr && gr->GetN() > 0) {
                            int n_points = gr->GetN();
                            double* x_vals = gr->GetX();

                            for (int i = 0; i < n_points; ++i) {
                                double x_val = x_vals[i];
                                
                                // Check each individual point against the floor limit
                                if (x_val > floor_limit) {
                                    if (x_val < true_min_x) {
                                        true_min_x = x_val;
                                    }
                                    if (x_val > true_max_x) {
                                        true_max_x = x_val;
                                    }
                                }
                            }
                        } else if (gr) {
                            std::cout << "  Graph: \"" << gr->GetName() << "\" <-- WARNING: Graph is empty (0 points)!" << std::endl;
                        }
                    }

                    if (true_min_x < INT_MAX && true_max_x > INT_MIN && true_max_x >= true_min_x) {
                        double safety_buffer = (true_max_x > true_min_x) ? (true_max_x - true_min_x) * 0.05 : 50.0; 

                        dynamic_lower_bound = true_min_x - safety_buffer;
                        dynamic_upper_bound = true_max_x + safety_buffer;
                    } else {
                        std::cout << "  -> WARNING: No valid points found above " << floor_limit << ". Falling back to default limits." << std::endl;
                    }
                }

                xAxis->SetLimits(dynamic_lower_bound, dynamic_upper_bound);
                xAxis->SetRangeUser(dynamic_lower_bound, dynamic_upper_bound);
                xAxis->SetTitle(x_label.c_str());
            }

            TAxis* yAxis = mg->GetHistogram()->GetYaxis();
            if (yAxis) {
                double dynamic_ymin = 0.0;
                double dynamic_ymax = 1.0;

                if (mg->GetListOfGraphs()) {
                    double true_min_y = INT_MAX;
                    double true_max_y = INT_MIN;

                    TIter next(mg->GetListOfGraphs());
                    TObject* gr_obj = nullptr;
                    while ((gr_obj = next())) {
                        auto* gr = dynamic_cast<TGraph*>(gr_obj);
                        if (!gr || gr->GetN() <= 0) continue;

                        // Try to cast to error-bearing graph subclasses
                        auto* gr_err = dynamic_cast<TGraphErrors*>(gr);
                        auto* gr_asymm = dynamic_cast<TGraphAsymmErrors*>(gr);

                        int n_points = gr->GetN();
                        double* x_vals = gr->GetX();
                        double* y_vals = gr->GetY();

                        for (int i = 0; i < n_points; ++i) {
                            // Skip y-values for any high-voltage points below or equal to the floor
                            if (x_vals && x_vals[i] <= floor_limit) {
                                continue; 
                            }

                            double val_y = y_vals[i];
                            double low_y = val_y;
                            double high_y = val_y;

                            if (gr_asymm) {
                                low_y  -= gr_asymm->GetErrorYlow(i);
                                high_y += gr_asymm->GetErrorYhigh(i);
                            } else if (gr_err) {
                                double err_y = gr_err->GetErrorY(i);
                                low_y  -= err_y;
                                high_y += err_y;
                            }

                            if (low_y < true_min_y) {
                                true_min_y = low_y;
                            }
                            if (high_y > true_max_y) {
                                true_max_y = high_y;
                            }
                        }
                    }

                    if (true_min_y < INT_MAX && true_max_y > INT_MIN && true_max_y > true_min_y) {
                        double safety_buffer = (true_max_y - true_min_y) * 0.05;

                        dynamic_ymin = std::max(true_min_y - safety_buffer, 0.0);
                        dynamic_ymax = true_max_y + safety_buffer;
                    }
                }

                yAxis->SetLimits(dynamic_ymin, dynamic_ymax);
                yAxis->SetRangeUser(dynamic_ymin, dynamic_ymax);
                yAxis->SetTitle(y_label.c_str());
            }
        }

        if (auto named_obj = dynamic_cast<TNamed*>(obj)) {
            named_obj->SetTitle(title.c_str());
        }

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        float start_x = 0.21;
        float start_y = 0.85;
        float offset_y = 0.04;
        drawATLASLabel(start_x, start_y, "Work in Progress");
        bool title_drawn = drawPlotTitle(obj, start_x, start_y - offset_y);

        drawATLASLegend(obj, start_x, title_drawn ? start_y - 2 * offset_y : start_y - offset_y, 13);
    }

    void styleStripDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(0.05);
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.14);

        applyATLASStyle(obj, canvas);

        auto h1 = dynamic_cast<TH1*>(obj);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(2.0);
        h1->SetLineStyle(1);

        // Format: TColor::GetColorTransparent(Color_Index, Alpha_Opacity_From_0_to_1)
        Int_t light_blue_transparent = TColor::GetColorTransparent(kAzure + 7, 0.30);
        h1->SetFillColor(light_blue_transparent);
        h1->SetFillStyle(1001); // 1001 = Solid fill style

        h1->Draw("HIST");

        canvas->Modified();
        canvas->Update();

        drawATLASLabel(0.21, 0.86, "Work in Progress");
        drawPlotTitle(obj, 0.21, 0.82);
    }

    void styleToTDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(0.05);
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.14);

        applyATLASStyle(obj, canvas);

        auto h1 = dynamic_cast<TH1*>(obj);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(2.0);
        h1->SetLineStyle(1);

        Int_t light_green_transparent = TColor::GetColorTransparent(kGreen + 2, 0.30);
        h1->SetFillColor(light_green_transparent);
        h1->SetFillStyle(1001);

        h1->Draw("HIST");

        canvas->Modified();
        canvas->Update();

        drawATLASLabel(0.90, 0.90, "Work in Progress", 33);
        drawPlotTitle(obj, 0.90, 0.86, 33);
    }

    void styleToTCombinedDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);
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

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            0.92, 0.82,               // Coordinates for the header box
            "Work in Progress",       // Status string
            plot_title,               // Title string
            32,                       // Alignment
            kWhite, 0.70f,            // 85% semi-transparent white background
            kBlack, 1,                // Black 1px border line
            0.01                      // Inner padding
        );

        double legend_y = header ? header->GetY1NDC() - 0.02 : 0.70;
        drawATLASLegend(obj, 0.92, legend_y, 33);

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
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(cl->InheritsFrom(TH2::Class()) ? 0.16 : 0.05);
        canvas->SetTopMargin(0.06);
        canvas->SetBottomMargin(0.14);

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
        
        canvas->Modified();
        canvas->Update();

        std::string plot_title = obj ? obj->GetTitle() : "";
        TPaveText* header = drawATLASHeaderBlock(
            0.21, 0.82,               // Coordinates for the header box
            "Work in Progress",       // Status string
            plot_title,               // Title string
            12,                       // Alignment
            kWhite, 0.70f,            // 85% semi-transparent white background
            kBlack, 1,                // Black 1px border line
            0.01                      // Inner padding
        );

        double legend_y = header ? header->GetY1NDC() - 0.02 : 0.70;
        drawATLASLegend(obj, 0.21, legend_y, 33);
    }

    static const std::vector<std::pair<PlotCategory, StylerFnPtr>> styler_map = {
        {PlotCategory::EfficiencyVsHV,          &styleEfficiencyVsHV},
        {PlotCategory::MeanClusterSizeVsHV,     &styleEfficiencyVsHV},
        {PlotCategory::NoiseRateVsHV,           &styleEfficiencyVsHV},
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