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
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>
#include <TVirtualPad.h>
#include <TMath.h>

#include "plotting/plotStyler.hpp"

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
            kAzure + 1,    // Sharp Blue
            kRed + 1,      // Deep Red
            kGreen + 2,    // Dark Green
            kOrange + 7    // Burnt Orange
        };

        void drawATLASLabel(double x_left, double y_top, const std::string& status) {
            double left_margin = 0.16;
            double top_margin = 0.06;
            if (gPad) {
                left_margin = gPad->GetLeftMargin();
                top_margin = gPad->GetTopMargin();
            }

            // Anchor exactly inside the active axis frame workspace
            double final_x = left_margin + x_left;
            double final_y = (1.0 - top_margin) - y_top;

            TLatex l;
            l.SetNDC();
            l.SetTextFont(72);  // Bold-Italic Helvetica for "ATLAS"
            l.SetTextColor(kBlack);
            l.SetTextSize(0.04);
            l.DrawLatex(final_x, final_y, "ATLAS");

            if (!status.empty()) {
                TLatex s;
                s.SetNDC();
                s.SetTextFont(42);  // Regular Helvetica
                s.SetTextColor(kBlack);
                s.SetTextSize(0.04);

                // 0.13 NDC offset ensures "ATLAS" and Status don't overlap
                s.DrawLatex(final_x + 0.10, final_y, status.c_str()); 
            }
        }

        bool drawPlotTitle(TObject* obj, double x_left, double y_top) {
            if (!obj || !gPad) {
                std::cout << "Error: Invalid object or pad for drawing plot title." << std::endl;
                return false;
            }

            std::string titleStr = obj->GetTitle();
            if (titleStr.empty()) return false;

            double left_margin = gPad->GetLeftMargin();
            double top_margin = gPad->GetTopMargin();

            double final_x = left_margin + x_left;
            double final_y = (1.0 - top_margin) - y_top;

            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.04);
            t.SetTextColor(kBlack);
            t.DrawLatex(final_x, final_y, titleStr.c_str());

            return true;
        }

        void drawATLASLegend(TObject* obj, double x_left, double y_top) {
        if (!obj || !gPad) return;

        TList* items_list = nullptr;
        bool is_graph = false;
        bool is_stack = false;

        // Check what container type we have
        if (auto* mg = dynamic_cast<TMultiGraph*>(obj)) {
            items_list = mg->GetListOfGraphs();
            is_graph = true;
        } else if (auto* stack = dynamic_cast<THStack*>(obj)) {
            items_list = stack->GetHists();
            is_stack = true;
        }

        if (!items_list || items_list->GetSize() == 0) return;

        double left_margin = gPad->GetLeftMargin();
        double top_margin = gPad->GetTopMargin();

        double final_x_left = left_margin + x_left;
        double final_y_top = (1.0 - top_margin) - y_top;

        double entry_height = 0.04; 
        double total_height = items_list->GetSize() * entry_height;
        double final_x_right = final_x_left + 0.25;
        double final_y_bottom = final_y_top - total_height;

        TLegend* leg = new TLegend(final_x_left, final_y_bottom, final_x_right, final_y_top);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);

        TIter next_item(items_list);
        TObject* child = nullptr;

        while ((child = next_item())) {
            std::string name = child->GetName();
            std::string label = child->GetTitle(); // Use the object's configured title by default

            // Fallback parsers for cleaner labels if titles aren't pretty
            if (is_graph) {
                if (name.find("layer0") != std::string::npos) label = "Layer 0";
                else if (name.find("layer1") != std::string::npos) label = "Layer 1";
                else if (name.find("layer2") != std::string::npos) label = "Layer 2";

                // Graphs get Points + Error Bars marker style
                leg->AddEntry(child, label.c_str(), "pe");
            } 
            else if (is_stack) {
                if (name.find("tot_eta1") != std::string::npos) label = "#eta_{1} Side";
                else if (name.find("tot_eta2") != std::string::npos) label = "#eta_{2} Side";
                
                // Filled histograms get Line + Fill box style
                leg->AddEntry(child, label.c_str(), "lf");
            }
        }

        leg->Draw();
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

            if (gPad) {
                gPad->SetTickx(1); // 1 = Draw ticks on top side
                gPad->SetTicky(1); // 1 = Draw ticks on right side
            }

            // Safely extract axes depending on the object type
            if (auto h2 = dynamic_cast<TH2*>(obj)) {
                xAxis = h2->GetXaxis();
                yAxis = h2->GetYaxis();
                zAxis = h2->GetZaxis();
                h2->SetStats(0);
                adjustDynamicCB(h2, pad); // Adjust color bar dynamically based on max value
            } else if (auto h = dynamic_cast<TH1*>(obj)) {
                xAxis = h->GetXaxis();
                yAxis = h->GetYaxis();
                h->SetMarkerStyle(20);
                h->SetMarkerSize(1.0);
                h->SetLineColor(kBlack);
                h->SetLineWidth(2);
                h->SetStats(0);
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
            TAxis* xAxis = mg->GetHistogram()->GetXaxis();
            if (xAxis) {
                double floor_limit = 4500.0;
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
                            double gr_min = TMath::MinElement(gr->GetN(), gr->GetX());
                            double gr_max = TMath::MaxElement(gr->GetN(), gr->GetX());

                            if (gr_min > floor_limit && gr_min < true_min_x) {
                                true_min_x = gr_min;
                            }
                            if (gr_max > gr_min && gr_max > true_max_x) {
                                true_max_x = gr_max;
                            }
                        }
                    }

                    if (true_max_x < INT_MAX && true_min_x > INT_MIN) {
                        double safety_buffer = (true_max_x - true_min_x) * 0.05;

                        dynamic_lower_bound = true_min_x - safety_buffer;
                        dynamic_upper_bound = true_max_x + safety_buffer;
                    }
                }

                xAxis->SetLimits(dynamic_lower_bound, dynamic_upper_bound);
                xAxis->SetRangeUser(dynamic_lower_bound, dynamic_upper_bound);
                xAxis->SetTitle(x_label.c_str());
            }

            TAxis* yAxis = mg->GetHistogram()->GetYaxis();
            if (yAxis) {
                double floor_limit = 0.0;
                double dynamic_ymin = floor_limit;
                double dynamic_ymax = floor_limit;

                if (mg->GetListOfGraphs()) {
                    double true_min_y = INT_MAX;
                    double true_max_y = INT_MIN;

                    TIter next(mg->GetListOfGraphs());
                    TObject* gr_obj = nullptr;
                    while ((gr_obj = next())) {
                        auto* gr = dynamic_cast<TGraph*>(gr_obj);
                        if (gr && gr->GetN() > 0) {
                            double gr_min = TMath::MinElement(gr->GetN(), gr->GetY());
                            double gr_max = TMath::MaxElement(gr->GetN(), gr->GetY());

                            if (gr_min > floor_limit && gr_min < true_min_y) {
                                true_min_y = gr_min;
                            }
                            if (gr_max > gr_min && gr_max > true_max_y) {
                                true_max_y = gr_max;
                            }
                        }
                    }

                    if (true_max_y < INT_MAX && true_min_y > INT_MIN) {
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

        drawATLASLabel(0.05, 0.07, "Work in Progress");
        bool title_drawn = drawPlotTitle(obj, 0.05, 0.12);

        drawATLASLegend(obj, 0.05, title_drawn ? 0.17 : 0.12);
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

        drawATLASLabel(0.05, 0.07, "Work in Progress");
        drawPlotTitle(obj, 0.05, 0.12);
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

        drawATLASLabel(0.46, 0.07, "Work in Progress");
        drawPlotTitle(obj, 0.46, 0.12);
    }

    void styleToTCombinedDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(0.05);
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.14);

        applyATLASStyle(obj, canvas);

        auto stack = dynamic_cast<THStack*>(obj);
        if (!stack) return;

        TList* hist_list = stack->GetHists();
        if (hist_list) {
            TIter next(hist_list);
            TH1* hist = nullptr;
            int index = 0;

            while ((hist = static_cast<TH1*>(next()))) {
                Color_t base_color = (index == 0) ? kBlue : kRed;

                // Sharp solid outline
                hist->SetLineColor(base_color);
                hist->SetLineWidth(2);
                hist->SetLineStyle(1);

                // Matching partially transparent fill (30% opacity)
                Int_t trans_color = TColor::GetColorTransparent(base_color, 0.30);
                hist->SetFillColor(trans_color);
                hist->SetFillStyle(1001);

                index++;
            }
        }
        stack->Draw("nostack hist");

        canvas->Modified();
        canvas->Update();

        drawATLASLabel(0.46, 0.07, "Work in Progress");
        bool title_drawn = drawPlotTitle(obj, 0.46, 0.12);

        // Draw the newly adjusted general legend
        drawATLASLegend(obj, 0.46, title_drawn ? 0.17 : 0.12);
    }

    void styleDefaultPlot(TObject* obj, TCanvas* canvas, TClass* cl) {
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(cl->InheritsFrom(TH2::Class()) ? 0.14 : 0.05);
        canvas->SetTopMargin(0.06);
        canvas->SetBottomMargin(0.14);

        if (auto gr = dynamic_cast<TGraph*>(obj)) {
            if (gr->GetHistogram()) {
                gr->GetHistogram()->SetStats(0);
            }
        }

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

        drawATLASLabel(0.05, 0.07, "Work in Progress");
        bool title_drawn = drawPlotTitle(obj, 0.05, 0.12);

        if (auto mg = dynamic_cast<TMultiGraph*>(obj)) {
            drawATLASLegend(mg, 0.05, title_drawn ? 0.17 : 0.12);
        }
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