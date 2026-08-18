#pragma once

#include <vector>
#include <string>
#include <Rtypes.h>

class TObject;
class TPad;
class TLatex;
class TPaveText;
class TLegend;

namespace PlotterHelpers {
namespace ATLASStyler {

    extern const std::vector<Color_t> ATLAS_PALETTE;

    std::vector<TLatex*> drawATLASLabel(float ndc_x, float ndc_y, const std::string& status, short alignment = 11);

    TLatex* drawPlotTitle(TObject* obj, float ndc_x, float ndc_y, short alignment = 11);

    TPaveText* drawATLASHeaderBlock(
        double ndc_x, double ndc_y,
        const std::string& status = "", const std::string& title = "", short alignment = 33,
        Color_t fillColor = 0, double fillAlpha = 0.70,
        Color_t borderColor = 1, int borderWidth = 1,
        double innerPadding = 0.01);

    TLegend* drawATLASLegend(TObject* obj, const std::vector<std::string>& legend_entries,
                             float ndc_x, float ndc_y, short alignment);

    void applyATLASStyle(TObject* obj, TPad* pad);

} // namespace ATLASStyler
} // namespace PlotterHelpers