#pragma once

#include <utility>

#include <TColor.h>

class TVirtualPad;
class TLatex;

namespace PlotterHelpers {
namespace Objects {

    void line(TVirtualPad* canvas, float x_start, float x_end, float y_start, float y_end,
              Color_t color = kBlack, int line_style = 1, int line_width = 1);

    void box(TVirtualPad* canvas, float x_start, float x_end, float y_start, float y_end,
             Color_t color = kWhite, float alpha = 1.0,
             int line_width = 1, int line_style = 1, Color_t line_color = kBlack, float line_alpha = 1.0);

    void hatchedRegion(TVirtualPad* canvas,
                       float x_start, float x_end, float y_start, float y_end,
                       int pattern = 3245, int line_width = 1,
                       Color_t color = kAzure + 2, float alpha = 0.5);

    std::pair<float, float> getTextSizeNDC(TLatex* latex_obj);

} // namespace Objects
} // namespace PlotterHelpers