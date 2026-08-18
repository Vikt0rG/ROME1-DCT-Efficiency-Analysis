#include <iostream>

#include <TLine.h>
#include <TPave.h>
#include <TBox.h>
#include <TLatex.h>
#include <TVirtualPad.h>
#include <TColor.h>

#include "objects.hpp"

namespace PlotterHelpers {
namespace Objects {

    void line(TVirtualPad* canvas, float x_start, float x_end, float y_start, float y_end,
              Color_t color, int line_style, int line_width) {
        if (!canvas) return;
        canvas->cd();
        TLine* line = new TLine(x_start, y_start, x_end, y_end);
        line->SetLineColor(color);
        line->SetLineStyle(line_style);
        line->SetLineWidth(line_width);
        line->Draw("same");
    }

    void box(TVirtualPad* canvas, float x_start, float x_end, float y_start, float y_end,
             Color_t color, float alpha,
             int line_width, int line_style, Color_t line_color, float line_alpha) {
        if (!canvas) return;
        canvas->cd();

        if (x_start >= x_end || y_start >= y_end) {
            std::cerr << "Error: Invalid coordinates for box.\n";
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

    void hatchedRegion(TVirtualPad* canvas, float x_start, float x_end, float y_start, float y_end,
                       int pattern, int line_width, Color_t color, float alpha) {
        if (!canvas) return;
        canvas->cd();

        if (x_start >= x_end || y_start >= y_end) {
            std::cerr << "Error: Invalid coordinates for hatched region.\n";
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

} // namespace Objects
} // namespace PlotterHelpers