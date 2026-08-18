#include <algorithm>
#include <cmath>

#include <TColor.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <THStack.h>
#include <TH2.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TPad.h>

#include "ATLASStyler.hpp"
#include "plotTypes.hpp"

namespace PlotterHelpers {
namespace ATLASStyler {

    const std::vector<Color_t> ATLAS_PALETTE = {
        static_cast<Color_t>(TColor::GetColor("#144d92")),
        static_cast<Color_t>(TColor::GetColor("#CF4446")),
        static_cast<Color_t>(TColor::GetColor("#1a8f3f")),
        static_cast<Color_t>(TColor::GetColor("#e28843"))
    };

    std::vector<TLatex*> drawATLASLabel(float ndc_x, float ndc_y,
        const std::string& status, short alignment) {
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

    TLatex* drawPlotTitle(TObject* obj, float ndc_x, float ndc_y, short alignment) {
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
        const std::string& status,
        const std::string& title,
        short alignment,
        Color_t fillColor, double fillAlpha,
        Color_t borderColor, int borderWidth,
        double innerPadding)
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

    // =========================================================================
    // ANONYMOUS NAMESPACE FOR PRIVATE HELPERS
    // =========================================================================
    namespace {
        struct SpacingSettings { double margin; double titleOffset; };

        int getLabelWidthDigits(double minVal, double maxVal, int nDivisions) {
            if (minVal == maxVal) return 1;

            // Characters needed for the integer part
            double maxAbs = std::max(std::abs(minVal), std::abs(maxVal));
            int intDigits = (maxAbs > 0) ? std::max(1, static_cast<int>(std::floor(std::log10(maxAbs))) + 1) : 1;

            // Number of primary ticks (ROOT default is Ndivisions = 510)
            int nPrimary = std::max(5, std::abs(nDivisions) % 100);
            double step = std::abs(maxVal - minVal) / nPrimary;

            // Characters needed for the decimal part
            int decDigits = 0;
            if (step > 0 && step < 1.0) {
                decDigits = std::clamp(static_cast<int>(std::ceil(-std::log10(step))), 0, 6);
            }

            return intDigits + decDigits;
        }

        SpacingSettings calculateAxisSpacing(double minVal, double maxVal, int nDivisions, AxisType type) {
            int digits = getLabelWidthDigits(minVal, maxVal, nDivisions);

            if (type == AxisType::X) {
                return {0.14, 1.10}; // X-Axis (Bottom margin & X-title offset)
            }

            if (type == AxisType::Y) { // Y-Axis (Left margin & Y-title offset)
                if (digits > 5)  return {0.18, 1.40};
                if (digits >= 4) return {0.16, 1.25};
                if (digits == 3) return {0.14, 1.10};
                if (digits == 2) return {0.12, 0.95};
                return {0.10, 0.80};
            }

            if (type == AxisType::Z) { // Colorbar (Right margin & Z-title offset)
                if (digits > 5)  return {0.18, 1.45};
                if (digits >= 4) return {0.17, 1.15};
                if (digits == 3) return {0.16, 1.00};
                if (digits == 2) return {0.15, 0.85};
                return {0.14, 0.70};
            }

            return {0.10, 1.00}; // Fallback
        }

        // Helper to isolate casting logic
        std::tuple<TAxis*, double, double> extractAxisData(TObject* obj, AxisType type) {
            if (!obj) return {nullptr, 0.0, 0.0};

            if (auto h2 = dynamic_cast<TH2*>(obj)) {
                if (type == AxisType::X && h2->GetXaxis()) return {h2->GetXaxis(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax()};
                if (type == AxisType::Y && h2->GetYaxis()) return {h2->GetYaxis(), h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax()};
                if (type == AxisType::Z && h2->GetZaxis()) return {h2->GetZaxis(), h2->GetMinimum(), h2->GetMaximum()};
            }
            else if (auto h1 = dynamic_cast<TH1*>(obj)) {
                if (type == AxisType::X && h1->GetXaxis()) return {h1->GetXaxis(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax()};
                if (type == AxisType::Y && h1->GetYaxis()) return {h1->GetYaxis(), h1->GetMinimum(), h1->GetMaximum()};
            }
            else if (auto stack = dynamic_cast<THStack*>(obj)) {
                if (type == AxisType::X && stack->GetXaxis()) return {stack->GetXaxis(), stack->GetXaxis()->GetXmin(), stack->GetXaxis()->GetXmax()};
                if (type == AxisType::Y && stack->GetYaxis()) return {stack->GetYaxis(), stack->GetMinimum(), stack->GetMaximum()};
            }
            else if (auto mg = dynamic_cast<TMultiGraph*>(obj)) {
                if (type == AxisType::X && mg->GetXaxis()) return {mg->GetXaxis(), mg->GetXaxis()->GetXmin(), mg->GetXaxis()->GetXmax()};
                if (type == AxisType::Y && mg->GetYaxis()) {
                    double minV = 1e9, maxV = -1e9;
                    if (mg->GetListOfGraphs()) {
                        for (TObject* gr_obj : *mg->GetListOfGraphs()) {
                            if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                                int n = gr->GetN();
                                if (n > 0) {
                                    minV = std::min(minV, *std::min_element(gr->GetY(), gr->GetY() + n));
                                    maxV = std::max(maxV, *std::max_element(gr->GetY(), gr->GetY() + n));
                                }
                            }
                        }
                    }
                    return {mg->GetYaxis(), minV, maxV};
                }
            }
            else if (auto g = dynamic_cast<TGraph*>(obj)) { // Handles TGraph, TGraphErrors, and TGraphAsymmErrors
                if (type == AxisType::X && g->GetXaxis()) return {g->GetXaxis(), g->GetXaxis()->GetXmin(), g->GetXaxis()->GetXmax()};
                if (type == AxisType::Y && g->GetYaxis()) {
                    int n = g->GetN();
                    if (n > 0) {
                        return {g->GetYaxis(), *std::min_element(g->GetY(), g->GetY() + n), *std::max_element(g->GetY(), g->GetY() + n)};
                    }
                }
            }
            return {nullptr, 0.0, 0.0};
        }

        // Master function to set pad margins and axis title offsets dynamically based on data range
        void setMargins(TObject* obj, TPad* pad, AxisType axisType) {
            if (!obj || !pad) return;

            // Safely extract the axis and true data limits
            auto [axis, minVal, maxVal] = extractAxisData(obj, axisType);
            if (!axis) return;

            // Default base margins
            pad->SetTopMargin(0.05);
            pad->SetRightMargin(0.05);

            // Calculate dynamic spacing
            SpacingSettings spacing = calculateAxisSpacing(minVal, maxVal, axis->GetNdivisions(), axisType);

            // Apply pad margins depending on axis orientation
            if (axisType == AxisType::Z)      pad->SetRightMargin(spacing.margin);
            else if (axisType == AxisType::Y) pad->SetLeftMargin(spacing.margin);
            else if (axisType == AxisType::X) pad->SetBottomMargin(spacing.margin);

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
    } // end anonymous namespace

    // =========================================================================

    void applyATLASStyle(TObject* obj, TPad* pad) {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetLineStyleString(11, "48 24");

        if (gPad) {
            gPad->SetTickx(1);
            gPad->SetTicky(1);
        }

        auto [xAxis, yAxis, zAxis] = extractAxis(obj, pad);

        if (xAxis) styleAxis(xAxis, pad, obj, AxisType::X, 0.04, 0.05, 510);
        if (yAxis) styleAxis(yAxis, pad, obj, AxisType::Y, 0.04, 0.05, 510);
        if (zAxis) styleAxis(zAxis, pad, nullptr, AxisType::Z, 0.04, 0.05, 510);

        if (auto h2 = dynamic_cast<TH2*>(obj)) adjustDynamicCB(h2, pad);
        if (zAxis) gStyle->SetPalette(kBird);

        pad->Modified();
        pad->Update();
    }

} // namespace ATLASStyler
} // namespace PlotterHelpers