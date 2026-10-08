#include <cmath>
#include <utility>
#include <optional>

#include <iostream>
#include <regex>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
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
#include <TFitResult.h>

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
        {"h2d_roi_map",             TH2::Class(),                  PlotCategory::RoI},
        {"h1d_cs_eta",              TH1::Class(),                  PlotCategory::CSDistribution},
        {"h1d_tot_eta",             TH1::Class(),                  PlotCategory::ToTDistribution},
        {"stack_tot_side_",         THStack::Class(),              PlotCategory::ToTCombSidesDistribution},
        {"stack_tot_layer",         THStack::Class(),              PlotCategory::ToTCombLayersDistribution},
        {"h1d_tof_layer",           TH1::Class(),                  PlotCategory::ToFDistribution},
        {"h2d_time_of_flight",      TH2::Class(),                  PlotCategory::ToFHeatmap},
        {"avg_time_of_flight",      TMultiGraph::Class(),          PlotCategory::AvgToFVsHV},
        {"time_resolution_layer",   TMultiGraph::Class(),          PlotCategory::TimeResolutionVsHV},
        {"time_resolution_strip",   TMultiGraph::Class(),          PlotCategory::TimeResolutionStripVsHV},
        {"track_eff",               TGraphAsymmErrors::Class(),    PlotCategory::Efficiency},
        {"eff",                     TGraphAsymmErrors::Class(),    PlotCategory::Efficiency},
        {"track_eff",               TMultiGraph::Class(),          PlotCategory::EfficiencyVsHV},
        {"eff",                     TMultiGraph::Class(),          PlotCategory::EfficiencyVsHV},
        {"avg_cluster_size",        TMultiGraph::Class(),          PlotCategory::MeanClusterSizeVsHV},
        {"rate_eta",                TMultiGraph::Class(),          PlotCategory::RateVsHV},
        {"rate_strips_eta",         TMultiGraph::Class(),          PlotCategory::RateStripVsHV},
        {"avg_tot_layer_eta",       TMultiGraph::Class(),          PlotCategory::AvgToTLayerVsHV},
        {"avg_tot_strip_eta",       TMultiGraph::Class(),          PlotCategory::AvgToTStripVsHV},
        {"avg_multiplicity",        TMultiGraph::Class(),          PlotCategory::AvgMultVsHV},
        {"stack_delay",             THStack::Class(),              PlotCategory::DelayDistribution},
        {"stack_avg_mult",          THStack::Class(),              PlotCategory::AvgMultStrip},
        {"stack_frac_mult",         THStack::Class(),              PlotCategory::FracMultStrip},
        {"h1d_delay",               TH1::Class(),                  PlotCategory::DelayDistribution},
        {"stack_avg_delay",         THStack::Class(),              PlotCategory::AvgDelayStrip},
        {"summary_v50_eff",         TGraphErrors::Class(),         PlotCategory::Efficiency},
        {"summary_v50_track_eff",   TGraphErrors::Class(),         PlotCategory::Efficiency},
        {"summary_vwp_eff",         TGraphErrors::Class(),         PlotCategory::Efficiency},
        {"summary_vwp_track_eff",   TGraphErrors::Class(),         PlotCategory::Efficiency}
    };

    static const std::vector<std::pair<PlotCategory, StylerFnPtr>> styler_map = {
        {PlotCategory::Efficiency,                  &styleEfficiency},
        {PlotCategory::EfficiencyVsHV,              &styleEfficiencyVsHV},
        {PlotCategory::MeanClusterSizeVsHV,         &styleAvgClusterSizeVsEff},
        {PlotCategory::RateVsHV,                    &styleRateVsEff},
        {PlotCategory::RateStripVsHV,               &styleRateStripVsHV},
        {PlotCategory::CSDistribution,              &styleCSDistribution},
        {PlotCategory::ToFDistribution,             &styleToFDistribution},
        {PlotCategory::ToFHeatmap,                  &styleToFHeatmap},
        {PlotCategory::AvgToFVsHV,                  &styleAvgToFVsHV},
        {PlotCategory::TimeResolutionVsHV,          &styleAvgTRVsHV},
        {PlotCategory::TimeResolutionStripVsHV,     &styleTRStripVsHV},
        {PlotCategory::AvgToTLayerVsHV,             &styleAvgClusterSizeVsEff},
        {PlotCategory::AvgToTStripVsHV,             &styleAvgToTVsHV},
        {PlotCategory::AvgMultVsHV,                 &styleAvgMulVsHV},
        {PlotCategory::StripDistribution,           &styleStripDistribution},
        {PlotCategory::StripDistributionCombined,   &styleStripDistributionCombined},
        {PlotCategory::RoI,                         &styleRoI},
        {PlotCategory::ToTDistribution,             &styleToTDistribution},
        {PlotCategory::ToTCombSidesDistribution,    &styleToTCombSidesDistribution},
        {PlotCategory::ToTCombLayersDistribution,   &styleToTCombLayersDistribution},
        {PlotCategory::AvgMultStrip,                &styleAvgMultStrip},
        {PlotCategory::FracMultStrip,               &styleFracMultStrip},
        {PlotCategory::AvgDelayStrip,               &styleAvgDelayStrip},
        {PlotCategory::DelayDistribution,           &styleDelayDistribution},
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

    std::optional<std::vector<std::string>> getLabelsFromGroupNames(const std::vector<std::string>& group_names) {
        std::vector<std::string> labels;
        for (const auto& group_name : group_names) {
            std::string label = group_name;

            if (label.rfind("group_", 0) != 0) return std::nullopt;

            std::string prefix = "group_";
            std::string clean_group = group_name.substr(prefix.length());

            std::smatch match;

            static const std::regex single_layer_re("layer[ _]?(\\d+)");
            static const std::regex lv_re("lv[ _]?(\\d+)");
            static const std::regex source_re("(?:source|filter)[ _]?(\\d+[._]\\d+|OFF)");
            static const std::regex mixture_re("mix(?:ture)?[ _]?([a-zA-Z0-9]+)");

            // A. Layer prefixes (e.g., "layer0" or "layer1")
            if (std::regex_search(clean_group, match, single_layer_re)) {
                labels.push_back("Layer " + match[1].str());
            }
            // B. Low voltage setting prefixes (e.g., "lv1", "lv2")
            else if (std::regex_search(clean_group, match, lv_re)) {
                labels.push_back("LV Setting " + match[1].str());
            }
            // C. Filter/Source prefixes
            else if (std::regex_search(clean_group, match, source_re)) {
                std::string extracted = match[1].str();
                if (extracted == "OFF") {
                    labels.push_back("Source OFF");
                } else {
                    std::replace(extracted.begin(), extracted.end(), '_', '.');
                    labels.push_back("Filter " + extracted);
                }
            }
            // D. Mixture prefixes (e.g., "mixSTD/mixECO1")
            else if (std::regex_search(clean_group, match, mixture_re)) {
                labels.push_back("Mixture " + match[1].str());
            }
            // Z. Default case: Use the cleaned group name as-is
            else {
                labels.push_back(clean_group);
            }

        }
        return labels;
    }

    std::string matchLabels(const std::string& name, const std::string& pattern) {
        std::regex re(pattern);
        std::smatch match;
        if (std::regex_search(name, match, re)) return match[1].str();
        return "";
    }

    std::tuple<std::vector<std::string>, std::string, std::string, std::vector<std::string>> compilePlotLabels(
        const std::string& metric_name,
        TObject* obj)
    {
        std::vector<std::string> out_title_lines;
        std::string out_xaxis = "High Voltage [V]", out_yaxis;
        std::vector<std::string> title_parts;
        std::vector<std::string> title_parts_extended;

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
            out_yaxis = "Time Resolution [ns]";
        } else if (metric_name.find("time_of_flight") != std::string::npos) {
            out_yaxis = metric_name.find("avg_") != std::string::npos ? "#LTToF#GT [ns]" : "Time of Flight [ns]";
        } else if (metric_name.find("strip") != std::string::npos) {
            out_xaxis = "Strip Number"; out_yaxis = "Hits";
        } else {
            out_xaxis = "X-Axis [a.u.]"; out_yaxis = "Value [a.u.]";
        }

        // -------------------------------------------------------------------------
        // Build Subtitle Context Pieces

        std::smatch match;

        // Heatmap Prefix
        if (metric_name.find("h2d_") == 0) {
            title_parts.push_back("Heatmap");
        }

        // Mixture Context
        static const std::regex mixture_re("mix(?:ture)?[ _]?([a-zA-Z0-9]+)");
        if (std::regex_search(metric_name, match, mixture_re)) {
            std::string mixture = match[1].str();
            title_parts.push_back(mixture + " Mixture");
        }

        // Layer Context
        static const std::regex single_layer_re("layer(\\d+)");
        if (std::regex_search(metric_name, match, single_layer_re)) {
            title_parts.push_back("Layer " + match[1].str());
        }

        // Side Context
        static const std::regex side_re("eta1|eta2|_or_|_and_");
        if (std::regex_search(metric_name, match, side_re)) {
            std::string m = match.str(0);
            if (m == "eta1")       title_parts.push_back("Side #eta_{1}");
            else if (m == "eta2")  title_parts.push_back("Side #eta_{2}");
            else if (m == "_or_")  title_parts.push_back("OR(#eta_{1}, #eta_{2})");
            else if (m == "_and_") title_parts.push_back("AND(#eta_{1}, #eta_{2})");
        }

        // Layer Pair Context
        static const std::regex layer_pair_re("layer_(\\d)_(\\d)");
        if (std::regex_search(metric_name, match, layer_pair_re)) {
            std::string new_title = Form("#it{t}_{Layer %s} #minus #it{t}_{Layer %s}",
                                         match[1].str().c_str(),
                                         match[2].str().c_str());
            title_parts.push_back(new_title);
        }

        // High Voltage Context
        static const std::regex hv_re("_hv_(\\d+)");
        if (std::regex_search(metric_name, match, hv_re)) {
            title_parts.push_back("HV: " + match[1].str() + " V");
        }

        // Track Reconstruction & Beam Context
        static const std::regex reco_re("^(beam_)?(track_)?(avg_tot|avg_multiplicity|eff)_");
        if (std::regex_search(metric_name, match, reco_re)) {

            // match[1] is the (beam_) group
            if (match[1].matched) {
                title_parts_extended.push_back("Beam Spot Region");
            }

            // match[2] is the (track_) group
            title_parts_extended.push_back(match[2].matched ? "After Track Reco" : "Before Track Reco");
        }

        // Trigger Context
        static const std::regex trigger_re("external|rpc");
        if (std::regex_search(metric_name, match, trigger_re)) {
            std::string m = match.str(0);
            if (m == "external") title_parts_extended.push_back("External Trigger");
            else if (m == "rpc") title_parts_extended.push_back("RPC Coincidence");
        }

        auto joinParts = [](const std::vector<std::string>& parts, const std::string& delim) {
            std::string result;
            for (size_t i = 0; i < parts.size(); ++i) {
                result += parts[i];
                if (i < parts.size() - 1) result += delim;
            }
            return result;
        };

        // Assemble Extended Title (e.g., "After Track Reco: External Trigger)
        if (!title_parts_extended.empty()) {
            out_title_lines.push_back(joinParts(title_parts_extended, ": "));
        }
        // Assemble Title (e.g., "Layers 0 & 1: Side #eta_{1}")
        if (!title_parts.empty()) {
            out_title_lines.push_back(joinParts(title_parts, ": "));
        }

        // -------------------------------------------------------------------------
        /// Build legend entries if object is a container (TMultiGraph or THStack)
        std::vector<std::string> legend_entries;
        if (TMultiGraph* mg = dynamic_cast<TMultiGraph*>(obj)) {
            TIter next(mg->GetListOfGraphs());
            TObject* obj;

            while ((obj = next())) {

                if (auto gr = dynamic_cast<TGraph*>(obj)) {
                    std::string gr_name = gr->GetTitle();
                    std::string layer = "", strip = "", mixture = "";

                    mixture = matchLabels(gr_name, "mix(?:ture)?[ _]?(\\w+)");
                    layer = matchLabels(gr_name, "layer(\\d+)");
                    strip = matchLabels(gr_name, "strip(\\d+)");
                    std::string legend_entry;

                    if (!mixture.empty()) {
                        legend_entry = mixture + " Mixture";
                    }
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
            } else {

                static const std::regex layer_re("layer(\\d+)");
                auto words_begin = std::sregex_iterator(metric_name.begin(), metric_name.end(), layer_re);
                auto words_end = std::sregex_iterator();

                for (std::sregex_iterator i = words_begin; i != words_end; ++i) {
                    std::string layer = (*i).str(1);
                    legend_entries.push_back("Layer " + layer);
                }
            }
        }

        return std::make_tuple(out_title_lines, out_xaxis, out_yaxis, legend_entries);
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
            // DATA 2: Histograms (Handles TH1, TH2, TH3)
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

            if (axis_type == AxisType::Y) {
                if (auto st = dynamic_cast<THStack*>(obj)) {
                    st->SetMinimum(dynamic_min);
                    st->SetMaximum(dynamic_max);
                } else if (auto h1 = dynamic_cast<TH1*>(obj)) {
                    h1->SetMinimum(dynamic_min);
                    h1->SetMaximum(dynamic_max);
                } else if (auto mg = dynamic_cast<TMultiGraph*>(obj)) {
                    mg->SetMinimum(dynamic_min);
                    mg->SetMaximum(dynamic_max);
                }
            }
        } else {
            if (default_min.has_value() && default_max.has_value()) {
                axis->SetLimits(default_min.value(), default_max.value());
                axis->SetRangeUser(default_min.value(), default_max.value());
                if (axis_type == AxisType::Y) {
                    if (auto st = dynamic_cast<THStack*>(obj)) {
                        st->SetMinimum(default_min.value());
                        st->SetMaximum(default_max.value());
                    } else if (auto h1 = dynamic_cast<TH1*>(obj)) {
                        h1->SetMinimum(default_min.value());
                        h1->SetMaximum(default_max.value());
                    } else if (auto mg = dynamic_cast<TMultiGraph*>(obj)) {
                        mg->SetMinimum(default_min.value());
                        mg->SetMaximum(default_max.value());
                    }
                }
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

    TGraph* findEfficiencyGraphForObject(TObject* obj) {
        if (!obj) return nullptr;
        TFile* file = gDirectory ? gDirectory->GetFile() : nullptr;
        if (!file) return nullptr;

        std::string obj_name = obj->GetName();
        std::string obj_title = "";
        if (auto named = dynamic_cast<TNamed*>(obj)) {
            obj_title = named->GetTitle();
        }
        std::string full_id = obj_name + " " + obj_title;

        // Extract group token cleanly (e.g., "group_mixECO1" from "rate_strips_eta2_group_mixECO1_layer2")
        std::string target_group = "";
        std::smatch match;

        if (std::regex_search(full_id, match, std::regex("(group_[A-Za-z0-9_]+?)(?=_layer|$)"))) {
            target_group = match[1].str();
        } else {
            return nullptr;
        }

        auto getEffFromDir = [](TDirectory* group_dir) -> TGraph* {
            if (!group_dir) return nullptr;
            TDirectory* eff_dir = group_dir->GetDirectory("efficiency_analysis");
            if (!eff_dir) return nullptr;

            TObject* eff_obj = eff_dir->Get("eff_or_rpc");
            if (!eff_obj) return nullptr;

            if (auto g = dynamic_cast<TGraph*>(eff_obj)) return g;
            if (auto mg = dynamic_cast<TMultiGraph*>(eff_obj)) {
                if (mg->GetListOfGraphs() && mg->GetListOfGraphs()->GetSize() > 0) {
                    return dynamic_cast<TGraph*>(mg->GetListOfGraphs()->At(0));
                }
            }
            return nullptr;
        };

        TIter next_top(file->GetListOfKeys());
        TKey* top_key = nullptr;
        while ((top_key = static_cast<TKey*>(next_top()))) {
            TClass* cl = TClass::GetClass(top_key->GetClassName());
            if (!cl || !cl->InheritsFrom(TDirectory::Class())) continue;

            TDirectory* config_dir = dynamic_cast<TDirectory*>(top_key->ReadObj());
            if (!config_dir) continue;

            TDirectory* group_dir = config_dir->GetDirectory(target_group.c_str());
            if (group_dir) {
                TGraph* g = getEffFromDir(group_dir);
                if (g) return g;
            }
        }
        return nullptr;
    }

    void convertMultiGraphToEfficiencyScale(TMultiGraph* mg) {
        if (!mg || !mg->GetListOfGraphs()) return;

        TFile* file = gDirectory ? gDirectory->GetFile() : nullptr;
        if (!file) return;

        auto getEffFromDir = [](TDirectory* group_dir) -> TGraph* {
            if (!group_dir) return nullptr;
            TDirectory* eff_dir = group_dir->GetDirectory("efficiency_analysis");
            if (!eff_dir) return nullptr;

            TObject* eff_obj = eff_dir->Get("eff_or_rpc");
            if (!eff_obj) return nullptr;

            if (auto g = dynamic_cast<TGraph*>(eff_obj)) return g;
            if (auto eff_mg = dynamic_cast<TMultiGraph*>(eff_obj)) {
                if (eff_mg->GetListOfGraphs() && eff_mg->GetListOfGraphs()->GetSize() > 0) {
                    return dynamic_cast<TGraph*>(eff_mg->GetListOfGraphs()->At(0));
                }
            }
            return nullptr;
        };

        TIter next_gr(mg->GetListOfGraphs());
        TObject* gr_obj = nullptr;

        while ((gr_obj = next_gr())) {
            auto gr = dynamic_cast<TGraph*>(gr_obj);
            if (!gr) continue;

            std::string gr_name = gr->GetName();

            // Extract "group_<token>" at the very end of gr_name
            std::string target_group = "";
            std::smatch match;
            if (std::regex_search(gr_name, match, std::regex("(group_[A-Za-z0-9_]+)$"))) {
                target_group = match[1].str();
            }

            if (target_group.empty()) {
                std::cout << "[WARNING convertMultiGraph] Could not extract group token from end of name: " << gr_name << std::endl;
                continue;
            }

            TGraph* eff_graph = nullptr;
            TIter next_top(file->GetListOfKeys());
            TKey* top_key = nullptr;

            while ((top_key = static_cast<TKey*>(next_top()))) {
                TClass* cl = TClass::GetClass(top_key->GetClassName());
                if (!cl || !cl->InheritsFrom(TDirectory::Class())) continue;

                TDirectory* config_dir = dynamic_cast<TDirectory*>(top_key->ReadObj());
                if (!config_dir) continue;

                TDirectory* group_dir = config_dir->GetDirectory(target_group.c_str());
                if (group_dir) {
                    eff_graph = getEffFromDir(group_dir);
                    if (eff_graph) break;
                }
            }

            if (!eff_graph) {
                std::cout << "[WARNING convertMultiGraph] Could not find efficiency graph for group: " << target_group << std::endl;
                continue;
            }

            TF1* fit_func = nullptr;
            if (eff_graph->GetListOfFunctions() && !eff_graph->GetListOfFunctions()->IsEmpty()) {
                fit_func = dynamic_cast<TF1*>(eff_graph->GetListOfFunctions()->Last());
            }

            int n_points = gr->GetN();
            double* x_vals = gr->GetX();
            double* y_vals = gr->GetY();

            for (int i = 0; i < n_points; ++i) {
                double hv = x_vals[i];
                double eff = fit_func ? fit_func->Eval(hv) : eff_graph->Eval(hv);
                gr->SetPoint(i, eff * 100.0, y_vals[i]);
            }
        }
    }

    void addEfficiencyTopAxis(TGraph* eff_graph, TPad* pad, TMultiGraph* mg = nullptr) {
        if (!eff_graph || !pad) return;

        pad->cd();
        pad->Update();

        double x_min = pad->GetUxmin();
        double x_max = pad->GetUxmax();
        double y_max = pad->GetUymax();
        double y_min = pad->GetUymin();

        // Dynamically detect subdivisions from the main axis
        int n_subdivisions = 5;
        if (mg && mg->GetHistogram()) {
            TAxis* xAxis = mg->GetHistogram()->GetXaxis();
            if (xAxis) {
                int ndiv = xAxis->GetNdivisions();
                int n3 = ndiv / 10000;
                if (n3 > 0) n_subdivisions = n3;
            }
        }

        double* hv_vals = eff_graph->GetX();
        double* eff_vals = eff_graph->GetY();
        int n_points = eff_graph->GetN();

        TF1* fit_func = nullptr;
        if (eff_graph->GetListOfFunctions() && !eff_graph->GetListOfFunctions()->IsEmpty()) {
            fit_func = dynamic_cast<TF1*>(eff_graph->GetListOfFunctions()->Last());
        }

        TLatex* tex = new TLatex();
        tex->SetTextFont(42);

        double major_tick_len = (y_max - y_min) * 0.03;
        double minor_tick_len = (y_max - y_min) * 0.015;

        // Title at Top Right
        tex->SetTextAlign(31);
        tex->SetTextSize(0.035);
        tex->DrawLatex(x_max, y_max + major_tick_len * 2.5, "Efficiency OR(#eta_{1}, #eta_{2}) [%]");

        // Numbers (Centered above major ticks)
        tex->SetTextAlign(21);
        tex->SetTextSize(0.025);

        for (int i = 0; i < n_points; ++i) {
            double hv = hv_vals[i];
            double eff = fit_func ? fit_func->Eval(hv) : eff_vals[i];

            if (hv >= x_min && hv <= x_max) {
                TLine* tick = new TLine(hv, y_max, hv, y_max - major_tick_len);
                tick->SetLineColor(kBlack);
                tick->SetLineWidth(1);
                tick->Draw();

                std::string eff_str = Form("%.1f", eff * 100);
                double label_y = y_max + major_tick_len * 0.4;
                tex->DrawLatex(hv, label_y, eff_str.c_str());
            }

            if (i < n_points - 1) {
                double hv_next = hv_vals[i + 1];
                for (int sub = 1; sub < n_subdivisions; ++sub) {
                    double hv_sub = hv + (hv_next - hv) * (static_cast<double>(sub) / n_subdivisions);
                    if (hv_sub >= x_min && hv_sub <= x_max) {
                        TLine* minor_tick = new TLine(hv_sub, y_max, hv_sub, y_max - minor_tick_len);
                        minor_tick->SetLineColor(kBlack);
                        minor_tick->SetLineWidth(1);
                        minor_tick->Draw();
                    }
                }
            }
        }
    }

    void styleEfficiency(TObject* obj, TCanvas* canvas, TClass* cl) {

        constexpr double y_max_padding = 0.30;
        constexpr double y_min_padding = 0.05;

        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), obj);

        if (auto gr = dynamic_cast<TGraphErrors*>(obj)) {
            std::vector<std::string> group_names;
            for (int i{}; i < gr->GetN(); ++i) {
                group_names.push_back(gr->GetXaxis()->GetBinLabel(i + 1));
            }

            auto group_labels = getLabelsFromGroupNames(group_names);
            if (group_labels.has_value()) {
                std::vector<std::string> labels = group_labels.value();

                for (size_t i{}; i < labels.size(); ++i) {
                    gr->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
                }
            }
        }

        obj->Draw("AP0");

        applyATLASStyle(obj, canvas);

        if (auto gr = dynamic_cast<TGraphErrors*>(obj)) {
            gr->GetXaxis()->SetLabelFont(42);
            gr->GetXaxis()->SetLabelSize(0.06);

            double min_y = std::numeric_limits<double>::max();
            double max_y = std::numeric_limits<double>::lowest();

            for (int i = 0; i < gr->GetN(); ++i) {
                double y = gr->GetY()[i];
                double y_err = gr->GetEY()[i];

                if (y - y_err < min_y) min_y = y - y_err;
                if (y + y_err > max_y) max_y = y + y_err;
            }

            double y_range = max_y - min_y;
            if (y_range <= 0) y_range = 1.0;

            gr->SetMaximum(max_y + (y_range * y_max_padding));
            gr->SetMinimum(min_y - (y_range * y_min_padding));
            gr->SetMarkerStyle(55);
            gr->SetMarkerSize(2.0);
            gr->SetMarkerColor(kAzure + 2);
            gr->SetLineWidth(2);
            gr->SetLineColorAlpha(kAzure + 2, 0.7);
        }

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            title_lines,
            12,
            kWhite, 0.70,
            kBlack, 0,
            0.01
        );

        canvas->Modified();
        canvas->Update();
    }

    void styleEfficiencyVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        // Extract title and axis labels from the object's title string
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        // Extract side context for y-axis & remove it from the header block
        std::string side_suffix = "";
        std::vector<std::pair<std::string, std::string>> side_map = {
            {"Side #eta_{1}", "#eta_{1}"},
            {"Side #eta_{2}", "#eta_{2}"},
            {"OR(#eta_{1}, #eta_{2})", "OR(#eta_{1}, #eta_{2})"},
            {"AND(#eta_{1}, #eta_{2})", "AND(#eta_{1}, #eta_{2})"}
        };

        for (auto& line : title_lines) {
            for (const auto& [search_str, suffix] : side_map) {
                size_t pos = line.find(search_str);
                if (pos != std::string::npos) {
                    side_suffix = suffix;

                    // Remove the substring from the title line
                    line.erase(pos, search_str.length());

                    // Clean up any dangling delimiters (": ") left over from the erasure
                    while (line.find(":  :") != std::string::npos) line.replace(line.find(":  :"), 4, ": ");
                    while (line.find(": :") != std::string::npos) line.replace(line.find(": :"), 3, ":");
                    if (line.find(": ") == 0) line.erase(0, 2);
                    if (line.length() >= 2 && line.substr(line.length() - 2) == ": ") line.erase(line.length() - 2);
                    break;
                }
            }
        }

        title_lines.erase(std::remove(title_lines.begin(), title_lines.end(), ""), title_lines.end());

        // Override the Y-axis label with the extracted context
        if (!side_suffix.empty()) {
            y_label = "Efficiency " + side_suffix;
        }

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
        const std::vector<Color_t> palette = {
            kAzure + 2, kGreen + 2, kOrange + 10, kMagenta + 2, kYellow - 3, kCyan - 4
        };
        if (mg && mg->GetListOfGraphs()) {

            double global_min_x = 1e9;
            double global_max_x = -1e9;
            TIter next_bounds(mg->GetListOfGraphs());
            TObject* gr_obj_bounds;
            while ((gr_obj_bounds = next_bounds())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj_bounds)) {
                    if (gr->GetN() > 0) {
                        global_min_x = std::min(global_min_x, TMath::MinElement(gr->GetN(), gr->GetX()));
                        global_max_x = std::max(global_max_x, TMath::MaxElement(gr->GetN(), gr->GetX()));
                    }
                }
            }

            if (global_min_x > global_max_x) {
                global_min_x = 4500.0;
                global_max_x = 6000.0;
            }

            double x_min_band = std::max(4500.0, global_min_x - 50.0);
            double x_max_band = global_max_x + 50.0;

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
                    gr->SetFillColorAlpha(color, 0.25);
                    gr->SetFillStyle(1001);

                    color_idx++;

                    // Extract the fitted sigmoid function from the graph's list of functions
                    TF1* sigmoid = dynamic_cast<TF1*>(gr->GetListOfFunctions()->First());
                    if (!sigmoid) continue;

                    sigmoid->SetRange(x_min_band, x_max_band);

                    // Sigmoid fit band calculation
                    double p0 = sigmoid->GetParameter(0);
                    double p1 = sigmoid->GetParameter(1);
                    double p2 = sigmoid->GetParameter(2);

                    double ep0 = sigmoid->GetParError(0);
                    double ep1 = sigmoid->GetParError(1);
                    double ep2 = sigmoid->GetParError(2);

                    int n_band_points = 200;
                    double step = (x_max_band - x_min_band) / n_band_points;

                    TGraphErrors* fit_band = new TGraphErrors(n_band_points);
                    double sigma_multiplier = 1.0;

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

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            title_lines,              // Title string
            12,                       // Alignment
            kWhite, 0.70,             // semi-transparent white background
            kBlack, 1,                // Black 1px border line
            0.01                      // Inner padding
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
        drawATLASLegend(obj, legend_entries, 0.18, legend_y, 13, "pef");

        canvas->Modified();
        canvas->Update();
    }

    void styleAvgClusterSizeVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        // Extract title and axis labels from the object's title string
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

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
        const std::vector<Color_t> palette = {
            kAzure + 2,
            kGreen + 2,
            kOrange + 10,
            kMagenta + 2,
            kYellow - 3,
            kCyan - 4
        };
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

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            title_lines,              // Title string
            12,                       // Alignment
            kWhite, 0.70,             // semi-transparent white background
            kBlack, 0,                // No border line
            0.01                      // Inner padding
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
        drawATLASLegend(obj, legend_entries, ndc_x0 + 0.03, legend_y, 13);

        canvas->Modified();
        canvas->Update();
    }

    void styleAvgClusterSizeVsEff(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        if (!mg) return;

        convertMultiGraphToEfficiencyScale(mg);

        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        canvas->cd();
        mg->Draw("AP0Z");

        if (mg->GetHistogram()) {
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                xAxis->SetLimits(50.0, 100.0);
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 50.0, .x_max = 100.0});
                xAxis->SetTitle("Efficiency OR(#eta_{1}, #eta_{2}) [%]");
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 50.0, .x_max = 100.0});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        const std::vector<Color_t> palette = {
            kAzure + 2, kGreen + 2, kOrange + 10, kMagenta + 2, kYellow - 3, kCyan - 4
        };
        if (mg->GetListOfGraphs()) {
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

        applyATLASStyle(obj, canvas);
        canvas->SetTickx(1);
        canvas->SetTicky(1);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03,
            ndc_y0 - 0.09,
            "Work in Progress",
            title_lines,
            12,
            kWhite, 0.70,
            kBlack, 0,
            0.01
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
        drawATLASLegend(obj, legend_entries, ndc_x0 + 0.03, legend_y, 13);

        canvas->Modified();
        canvas->Update();
    }

    void styleRateVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        // Extract title and axis labels from the object's title string
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        obj->Draw("AP0Z");

        // Set axis ranges and labels
        if (mg && mg->GetHistogram()) {
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 4550.0});
                xAxis->SetTitle(x_label.c_str());
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 4550.0});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        // Set color and marker style for the graphs in the multigraph
        const std::vector<Color_t> palette = {
            kAzure + 2,
            kGreen + 2,
            kOrange + 10,
            kMagenta + 2,
            kYellow - 3,
            kCyan - 4
        };
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

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            title_lines,              // Title string
            12,                       // Alignment
            kWhite, 0.70,             // semi-transparent white background
            kBlack, 0,                // No border line
            0.01                      // Inner padding
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
        drawATLASLegend(obj, legend_entries, ndc_x0 + 0.03, legend_y, 13);

        canvas->Modified();
        canvas->Update();
    }

    void styleRateVsEff(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        if (!mg) return;

        convertMultiGraphToEfficiencyScale(mg);

        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        obj->Draw("AP0Z");

        // Set axis ranges and labels
        if (mg->GetHistogram()) {
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 0.0, .x_max = 100.0});
                xAxis->SetTitle(x_label.c_str());
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 0.0, .x_max = 100.0});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        // Set color and marker style for the graphs in the multigraph
        const std::vector<Color_t> palette = {
            kAzure + 2, kGreen + 2, kOrange + 10, kMagenta + 2, kYellow - 3, kCyan - 4
        };
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

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            title_lines,              // Title string
            12,                       // Alignment
            kWhite, 0.70,             // semi-transparent white background
            kBlack, 0,                // No border line
            0.01                      // Inner padding
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header ? header->GetY1NDC() - 0.04 : 0.70;
        drawATLASLegend(obj, legend_entries, ndc_x0 + 0.03, legend_y, 13);

        canvas->Modified();
        canvas->Update();
    }

    void styleRateStripVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        int n_graphs = (mg && mg->GetListOfGraphs()) ? mg->GetListOfGraphs()->GetSize() : 0;
        bool is_strip_plot = (n_graphs > 3);

        if (is_strip_plot) gStyle->SetPalette(kViridis);

        obj->Draw("APZ");

        if (mg && mg->GetHistogram()) {
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 4550.0});
                xAxis->SetTitle(x_label.c_str());
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 4550.0});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        applyATLASStyle(obj, canvas);
        canvas->SetTopMargin(0.12);

        if (is_strip_plot && n_graphs > 0) {
            canvas->SetRightMargin(0.16);

            // Find min and max strip indices present in this specific multigraph
            int min_strip = INT_MAX;
            int max_strip = INT_MIN;
            TIter pass1(mg->GetListOfGraphs());
            TObject* gr_obj_p1;
            while ((gr_obj_p1 = pass1())) {
                std::smatch match;
                std::string gr_title = gr_obj_p1->GetTitle();
                if (std::regex_search(gr_title, match, std::regex("Strip (\\d+)"))) {
                    int s_idx = std::stoi(match[1].str());
                    if (s_idx < min_strip) min_strip = s_idx;
                    if (s_idx > max_strip) max_strip = s_idx;
                }
            }
            if (min_strip > max_strip) { min_strip = 0; max_strip = STRIPS_PER_LAYER - 1; }
            int n_strips_active = max_strip - min_strip + 1;
            int n_colors = TColor::GetNumberOfColors();

            TIter next(mg->GetListOfGraphs());
            TObject* gr_obj;
            while ((gr_obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                    int strip_idx = min_strip;
                    std::smatch match;
                    std::string gr_title = gr->GetTitle();

                    if (std::regex_search(gr_title, match, std::regex("Strip (\\d+)"))) {
                        strip_idx = std::stoi(match[1].str());
                    }
                    strip_idx = std::max(min_strip, std::min(strip_idx, max_strip));

                    int bin = strip_idx - min_strip;
                    int color_idx = TColor::GetColorPalette((bin * n_colors) / n_strips_active);

                    gr->SetMarkerColor(color_idx);
                    gr->SetMarkerStyle(70);
                    gr->SetMarkerSize(1.8);
                    gr->SetLineColor(color_idx);
                    gr->SetLineWidth(2.0);
                }
            }

            TH2D* dummy_z = new TH2D(Form("dummy_z_%p", mg), "", 1, -2000, -1000, 1, -2000, -1000);
            dummy_z->SetDirectory(nullptr);
            dummy_z->SetBinContent(1, 1, 0.0);
            dummy_z->SetMinimum(min_strip);
            dummy_z->SetMaximum(max_strip + 1);
            dummy_z->SetContour(n_strips_active);

            TAxis* zAxis = dummy_z->GetZaxis();
            zAxis->SetTitle("Strip Number");
            zAxis->SetTitleOffset(1.0);
            zAxis->SetTitleSize(0.05);
            zAxis->SetLabelSize(0.04);

            dummy_z->Draw("COL Z SAME");
        }

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            title_lines,
            12,
            kWhite, 0.70,
            kBlack, 0,
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

        canvas->SetTickx(0);

        canvas->RedrawAxis();
        canvas->Modified();
        canvas->Update();

        TGraph* eff_graph = findEfficiencyGraphForObject(obj);
        if (eff_graph) {
            addEfficiencyTopAxis(eff_graph, canvas, mg);
            canvas->Update();
        }
    }

    void styleAvgToFVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {
        // Extract title and axis labels from the object's title string
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

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
        const std::vector<Color_t> palette = {
            kAzure + 2,
            kGreen + 2,
            kOrange + 10,
            kMagenta + 2,
            kYellow - 3,
            kCyan - 4
        };
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

        applyATLASStyle(obj, canvas);

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.03,
            ndc_y0 - 0.09,            // Coordinates for the header box
            "Work in Progress",       // Status string
            title_lines,              // Title string
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

    void styleAvgTRVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        if (!mg) return;
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        canvas->cd();
        mg->Draw("AP0Z");
        canvas->Update();

        double fit_x_max = 6000.0;

        // Set axis ranges and labels
        if (TH1* frame = mg->GetHistogram()) {
            if (TAxis* xAxis = frame->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 5200.0});
                xAxis->SetTitle(x_label.c_str());
                fit_x_max = xAxis->GetXmax();
            }

            if (TAxis* yAxis = frame->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 5200.0});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        applyATLASStyle(obj, canvas);
        canvas->cd();

        const std::vector<Color_t> palette = {
            kAzure + 2, kGreen + 2, kOrange + 10, kMagenta + 2, kYellow - 3, kCyan - 4
        };

        if (mg->GetListOfGraphs()) {

            double global_max_x = 0.0;
            TIter next_max(mg->GetListOfGraphs());
            TObject* gr_obj_max;
            while ((gr_obj_max = next_max())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj_max)) {
                    global_max_x = std::max(global_max_x, TMath::MaxElement(gr->GetN(), gr->GetX()));
                }
            }

            double fit_x_min = 5000.0;
            double end_x = std::min(fit_x_max, global_max_x) + 50.0;

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

                    gr->SetFillColorAlpha(color, 0.25);
                    gr->SetFillStyle(1001);

                    // Prevent ROOT's global function registry from deleting fits during batch runs
                    std::string fit_name = Form("fit_pol1_%p_%d", (void*)gr, color_idx);
                    TF1* fit = new TF1(fit_name.c_str(), "pol1", fit_x_min, end_x);

                    TFitResultPtr r = gr->Fit(fit, "Q0SR", "", fit_x_min, end_x);

                    if (static_cast<int>(r) == 0) {
                        double p0 = r->Parameter(0);
                        double p1 = r->Parameter(1);
                        double ep0 = r->Error(0);
                        double ep1 = r->Error(1);
                        double cov01 = r->CovMatrix(0, 1);

                        // Account for excess scatter in the data
                        double chi2 = r->Chi2();
                        double ndf = r->Ndf();
                        double scatter_scale = (ndf > 0 && (chi2 / ndf) > 1.0) ? std::sqrt(chi2 / ndf) : 1.0;

                        int n_points = 200;
                        double step = (end_x - fit_x_min) / n_points;
                        TGraphErrors* band = new TGraphErrors(n_points);

                        for (int i = 0; i < n_points; ++i) {
                            double x = fit_x_min + i * step;
                            double y = p0 + p1 * x;

                            double var = (ep0 * ep0) + (x * x * ep1 * ep1) + (2.0 * x * cov01);

                            // Apply the 1-Sigma multiplier and the scatter scale
                            double err = scatter_scale * std::sqrt(std::max(0.0, var));

                            band->SetPoint(i, x, y);
                            band->SetPointError(i, 0, err);
                        }

                        band->SetFillColorAlpha(color, 0.25);
                        band->SetLineColor(color);
                        band->SetLineWidth(1);
                        band->Draw("E3");

                        fit->SetLineColor(color);
                        fit->SetLineStyle(2);
                        fit->SetLineWidth(2);
                        fit->Draw("SAME");
                    }

                    color_idx++;
                }
            }
            mg->Draw("P0Z");
        }

        canvas->RedrawAxis();
        canvas->Modified();
        canvas->Update();

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.03,
            ndc_y0 - 0.09,
            "Work in Progress",
            title_lines,
            32,
            kWhite, 0.70,
            kBlack, 0,
            0.01
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header->GetY1NDC() - 0.02;
        int alignment = 33;
        drawATLASLegend(obj, legend_entries, ndc_x0, legend_y, alignment, "pef");

        canvas->Modified();
        canvas->Update();
    }

    void styleAvgTRVsEff(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        if (!mg) return;

        convertMultiGraphToEfficiencyScale(mg);

        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        canvas->cd();
        mg->Draw("AP0Z");
        canvas->Update();

        double x_max_limit = 100.0; // The rightmost extrapolation boundary

        if (TH1* frame = mg->GetHistogram()) {
            if (TAxis* xAxis = frame->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 10.0, .x_max = 100.0});
                xAxis->SetTitle("Efficiency OR(#eta_{1}, #eta_{2}) [%]");
                x_max_limit = xAxis->GetXmax();
            }
            if (TAxis* yAxis = frame->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 10.0, .x_max = 100.0});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        applyATLASStyle(obj, canvas);
        canvas->SetTickx(1);
        canvas->SetTicky(1);
        canvas->cd();

        const std::vector<Color_t> palette = {
            kAzure + 2, kGreen + 2, kOrange + 10, kMagenta + 2, kYellow - 3, kCyan - 4
        };

        if (mg->GetListOfGraphs()) {
            double fit_x_min = 10.0;

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

                    gr->SetFillColorAlpha(color, 0.25);
                    gr->SetFillStyle(1001);

                    double graph_max_x = TMath::MaxElement(gr->GetN(), gr->GetX());
                    double fit_x_max = std::min(x_max_limit, graph_max_x);

                    std::string fit_name = Form("fit_pol1_eff_%p_%d", (void*)gr, color_idx);
                    TF1* fit = new TF1(fit_name.c_str(), "pol1", fit_x_min, x_max_limit);

                    TFitResultPtr r = gr->Fit(fit, "Q0SR", "", fit_x_min, fit_x_max);

                    if (static_cast<int>(r) == 0) {
                        double p0 = r->Parameter(0);
                        double p1 = r->Parameter(1);
                        double ep0 = r->Error(0);
                        double ep1 = r->Error(1);
                        double cov01 = r->CovMatrix(0, 1);

                        double chi2 = r->Chi2();
                        double ndf = r->Ndf();
                        double scatter_scale = (ndf > 0 && (chi2 / ndf) > 1.0) ? std::sqrt(chi2 / ndf) : 1.0;

                        int n_points = 200;
                        double step = (x_max_limit - fit_x_min) / n_points;
                        TGraphErrors* band = new TGraphErrors(n_points);

                        for (int i = 0; i < n_points; ++i) {
                            double x = fit_x_min + i * step;
                            double y = p0 + p1 * x;

                            double var = (ep0 * ep0) + (x * x * ep1 * ep1) + (2.0 * x * cov01);
                            double err = scatter_scale * std::sqrt(std::max(0.0, var));

                            band->SetPoint(i, x, y);
                            band->SetPointError(i, 0, err);
                        }

                        band->SetFillColorAlpha(color, 0.25);
                        band->SetLineColor(color);
                        band->SetLineWidth(1);
                        band->Draw("E3");

                        fit->SetLineColor(color);
                        fit->SetLineStyle(2);
                        fit->SetLineWidth(2);
                        fit->Draw("SAME");
                    }

                    color_idx++;
                }
            }
            mg->Draw("P0Z");
        }

        canvas->RedrawAxis();
        canvas->Modified();
        canvas->Update();

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.03,
            ndc_y0 - 0.09,
            "Work in Progress",
            title_lines,
            32,
            kWhite, 0.70,
            kBlack, 0,
            0.01
        );

        canvas->Modified();
        canvas->Update();

        double legend_y = header->GetY1NDC() - 0.02;
        int alignment = 33;
        drawATLASLegend(obj, legend_entries, ndc_x0, legend_y, alignment, "pef");

        canvas->Modified();
        canvas->Update();
    }

    void styleTRStripVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        int n_graphs = (mg && mg->GetListOfGraphs()) ? mg->GetListOfGraphs()->GetSize() : 0;
        bool is_strip_plot = (n_graphs > 3);

        if (is_strip_plot) gStyle->SetPalette(kViridis);

        obj->Draw("APZ");

        if (mg && mg->GetHistogram()) {
            if (TAxis* xAxis = mg->GetHistogram()->GetXaxis()) {
                setRange(mg, xAxis, AxisType::X, std::nullopt, std::nullopt, {.x_min = 5150.0, .y_max = 1.20});
                xAxis->SetTitle(x_label.c_str());
            }
            if (TAxis* yAxis = mg->GetHistogram()->GetYaxis()) {
                setRange(mg, yAxis, AxisType::Y, std::nullopt, std::nullopt, {.x_min = 5150.0, .y_max = 1.20});
                yAxis->SetTitle(y_label.c_str());
            }
        }

        applyATLASStyle(obj, canvas);

        if (is_strip_plot && n_graphs > 0) {
            canvas->SetTopMargin(0.12);
            canvas->SetRightMargin(0.16);

            // Find min and max strip indices present in this specific multigraph
            int min_strip = INT_MAX;
            int max_strip = INT_MIN;
            TIter pass1(mg->GetListOfGraphs());
            TObject* gr_obj_p1;
            while ((gr_obj_p1 = pass1())) {
                std::smatch match;
                std::string gr_title = gr_obj_p1->GetTitle();
                if (std::regex_search(gr_title, match, std::regex("Strip (\\d+)"))) {
                    int s_idx = std::stoi(match[1].str());
                    if (s_idx < min_strip) min_strip = s_idx;
                    if (s_idx > max_strip) max_strip = s_idx;
                }
            }
            if (min_strip > max_strip) { min_strip = 0; max_strip = STRIPS_PER_LAYER - 1; }
            int n_strips_active = max_strip - min_strip + 1;
            int n_colors = TColor::GetNumberOfColors();

            TIter next(mg->GetListOfGraphs());
            TObject* gr_obj;
            while ((gr_obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                    int strip_idx = min_strip;
                    std::smatch match;
                    std::string gr_title = gr->GetTitle();

                    if (std::regex_search(gr_title, match, std::regex("Strip (\\d+)"))) {
                        strip_idx = std::stoi(match[1].str());
                    }
                    strip_idx = std::max(min_strip, std::min(strip_idx, max_strip));

                    int bin = strip_idx - min_strip;
                    int color_idx = TColor::GetColorPalette((bin * n_colors) / n_strips_active);

                    gr->SetMarkerColor(color_idx);
                    gr->SetMarkerStyle(70);
                    gr->SetMarkerSize(1.8);
                    gr->SetLineColor(color_idx);
                    gr->SetLineWidth(2.0);
                }
            }

            TH2D* dummy_z = new TH2D(Form("dummy_z_%p", mg), "", 1, -2000, -1000, 1, -2000, -1000);
            dummy_z->SetDirectory(nullptr);
            dummy_z->SetBinContent(1, 1, 0.0);
            dummy_z->SetMinimum(min_strip);
            dummy_z->SetMaximum(max_strip + 1);
            dummy_z->SetContour(n_strips_active);

            TAxis* zAxis = dummy_z->GetZaxis();
            zAxis->SetTitle("Strip Number");
            zAxis->SetTitleOffset(1.0);
            zAxis->SetTitleSize(0.05);
            zAxis->SetLabelSize(0.04);

            dummy_z->Draw("COL Z SAME");
        }

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            title_lines,
            12,
            kWhite, 0.70,
            kBlack, 0,
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

        if (is_strip_plot) {
            canvas->SetTickx(0);
            canvas->RedrawAxis();
        }

        canvas->Modified();
        canvas->Update();

        if (is_strip_plot) {
            TGraph* eff_graph = findEfficiencyGraphForObject(obj);
            if (eff_graph) {
                addEfficiencyTopAxis(eff_graph, canvas, mg);
                canvas->Update();
            }
        }
    }

    void styleAvgToTVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        int n_graphs = (mg && mg->GetListOfGraphs()) ? mg->GetListOfGraphs()->GetSize() : 0;
        bool is_strip_plot = (n_graphs > 3);

        if (is_strip_plot) gStyle->SetPalette(kViridis);

        obj->Draw("APZ");

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

        applyATLASStyle(obj, canvas);

        if (is_strip_plot && n_graphs > 0) {
            canvas->SetTopMargin(0.12);
            canvas->SetRightMargin(0.16);

            // Find min and max strip indices present in this specific multigraph
            int min_strip = INT_MAX;
            int max_strip = INT_MIN;
            TIter pass1(mg->GetListOfGraphs());
            TObject* gr_obj_p1;
            while ((gr_obj_p1 = pass1())) {
                std::smatch match;
                std::string gr_title = gr_obj_p1->GetTitle();
                if (std::regex_search(gr_title, match, std::regex("Strip (\\d+)"))) {
                    int s_idx = std::stoi(match[1].str());
                    if (s_idx < min_strip) min_strip = s_idx;
                    if (s_idx > max_strip) max_strip = s_idx;
                }
            }
            if (min_strip > max_strip) { min_strip = 0; max_strip = STRIPS_PER_LAYER - 1; }
            int n_strips_active = max_strip - min_strip + 1;
            int n_colors = TColor::GetNumberOfColors();

            // Recolor the graphs based on their actual relative position in the subset
            TIter next(mg->GetListOfGraphs());
            TObject* gr_obj;
            while ((gr_obj = next())) {
                if (auto gr = dynamic_cast<TGraph*>(gr_obj)) {
                    int strip_idx = min_strip;
                    std::smatch match;
                    std::string gr_title = gr->GetTitle();

                    if (std::regex_search(gr_title, match, std::regex("Strip (\\d+)"))) {
                        strip_idx = std::stoi(match[1].str());
                    }
                    strip_idx = std::max(min_strip, std::min(strip_idx, max_strip));

                    // Calculate the color index based on the strip index
                    int bin = strip_idx - min_strip;
                    int color_idx = TColor::GetColorPalette((bin * n_colors) / n_strips_active);

                    gr->SetMarkerColor(color_idx);
                    gr->SetMarkerStyle(70);
                    gr->SetMarkerSize(1.8);
                    gr->SetLineColor(color_idx);
                    gr->SetLineWidth(2.0);
                }
            }

            // Draw dummy Z-axis matching the active strip bounds
            TH2D* dummy_z = new TH2D(Form("dummy_z_%p", mg), "", 1, -2000, -1000, 1, -2000, -1000);
            dummy_z->SetDirectory(nullptr);
            dummy_z->SetBinContent(1, 1, 0.0);
            dummy_z->SetMinimum(min_strip);
            dummy_z->SetMaximum(max_strip + 1);
            dummy_z->SetContour(n_strips_active);

            TAxis* zAxis = dummy_z->GetZaxis();
            zAxis->SetTitle("Strip Number");
            zAxis->SetTitleOffset(1.0);
            zAxis->SetTitleSize(0.05);
            zAxis->SetLabelSize(0.04);

            dummy_z->Draw("COL Z SAME");
        }

        canvas->Modified();
        canvas->Update();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            title_lines,
            12,
            kWhite, 0.70,
            kBlack, 0,
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

        if (is_strip_plot) {
            canvas->SetTickx(0);
            canvas->RedrawAxis();
        }

        canvas->Modified();
        canvas->Update();

        if (is_strip_plot) {
            TGraph* eff_graph = findEfficiencyGraphForObject(obj);
            if (eff_graph) {
                addEfficiencyTopAxis(eff_graph, canvas, mg);
                canvas->Update();
            }
        }
    }

    void styleAvgMulVsHV(TObject* obj, TCanvas* canvas, TClass* cl) {

        auto mg = dynamic_cast<TMultiGraph*>(obj);
        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), mg);

        // Differentiate between a 24-strip plot and a 3-layer plot
        bool is_strip_plot = (mg && mg->GetListOfGraphs() && mg->GetListOfGraphs()->GetSize() > 10);

        if (is_strip_plot) gStyle->SetPalette(kViridis);

        obj->Draw("APZ");

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

        applyATLASStyle(obj, canvas);

        // Add a strip colorbar
        if (is_strip_plot) {
            canvas->SetTopMargin(0.12);
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

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            title_lines,
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

        if (is_strip_plot) {
            canvas->SetTickx(0);
            canvas->RedrawAxis();
        }

        canvas->Modified();
        canvas->Update();

        if (is_strip_plot) {
            TGraph* eff_graph = findEfficiencyGraphForObject(obj);
            if (eff_graph) {
                addEfficiencyTopAxis(eff_graph, canvas, mg);
                canvas->Update();
            }
        }
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

        enforceIntegerMinorTicks(h1->GetXaxis());

        drawATLASLabel(0.21, 0.86, "Work in Progress");
        drawPlotTitle(obj, 0.21, 0.82);
    }

    void styleStripDistributionCombined(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);
        if (!stack) return;

        stack->Draw("nostack hist");

        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetName(), stack);

        TIter next(stack->GetHists());
        TH1* hist = nullptr;
        int index = 0;

        while ((hist = static_cast<TH1*>(next()))) {
            Color_t base_color = (index == 0) ? kRed + 1 : kAzure - 3;

            hist->SetLineColor(base_color);
            hist->SetLineWidth(2);
            hist->SetLineStyle(1);

            Int_t trans_color = TColor::GetColorTransparent(base_color, 0.30);
            hist->SetFillColor(trans_color);
            hist->SetFillStyle(1001);

            for (int i = 1; i < hist->GetNbinsX(); ++i) {
                double x = hist->GetXaxis()->GetBinUpEdge(i);

                double y_left = hist->GetBinContent(i);
                double y_right = hist->GetBinContent(i + 1);

                double y_max_line = std::min(y_left, y_right);

                if (y_max_line > 0) {
                    TLine* edge = new TLine(x, 0.0, x, y_max_line);
                    edge->SetLineColorAlpha(base_color, 0.55);
                    edge->SetLineWidth(2);
                    edge->Draw();
                }
            }

            index++;
        }

        if (stack->GetXaxis()) stack->GetXaxis()->SetTitle(x_label.c_str());
        if (stack->GetYaxis()) stack->GetYaxis()->SetTitle(y_label.c_str());

        applyATLASStyle(obj, canvas);
        enforceIntegerMinorTicks(stack->GetXaxis());

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.09,
            "Work in Progress",
            title_lines,
            12,
            kWhite, 0.70,
            kBlack, 0,
            0.01
        );

        drawATLASLegend(obj, legend_entries, ndc_x0 + 0.03, ndc_y0 - 0.2, 13);

        canvas->Modified();
        canvas->Update();
    }

    void styleCSDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto h1 = dynamic_cast<TH1*>(obj);
        if (!h1) return;

        double x_max = h1->GetXaxis()->GetXmax();
        double y_max = h1->GetMaximum() * 1.15;

        TH1F* frame = canvas->DrawFrame(0.0, 0.0, x_max, y_max);
        frame->SetTitle(h1->GetTitle());
        frame->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        frame->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());

        applyATLASStyle(frame, canvas);

        h1->SetFillColorAlpha(kOrange - 2, 0.75);
        h1->SetFillStyle(1001);
        h1->SetLineColor(kOrange + 7);
        h1->SetLineWidth(2);
        h1->SetLineStyle(1);

        h1->Draw("HIST SAME");

        for (int i = 1; i < h1->GetNbinsX(); ++i) {
            double x = h1->GetXaxis()->GetBinUpEdge(i);

            double y_left = h1->GetBinContent(i);
            double y_right = h1->GetBinContent(i + 1);

            double y_max_line = std::min(y_left, y_right);

            if (y_max_line > 0) {
                TLine* edge = new TLine(x, 0.0, x, y_max_line);
                edge->SetLineColorAlpha(kOrange + 7, 0.55);
                edge->SetLineWidth(1);
                edge->Draw();
            }
        }

        TF1* fit = h1->GetFunction("cs_poisson");

        if (!fit) {
            std::cerr << "No Poisson fit found for histogram: " << h1->GetName() << std::endl;
            return;
        }

        double p1 = fit->GetParameter(1);
        double ep1 = fit->GetParError(1);

        canvas->RedrawAxis();

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TLegend* fit_legend = new TLegend(ndc_x0 - 0.40, ndc_y0 - 0.3, ndc_x0 - 0.03, ndc_y0 - 0.13);
        fit_legend->SetTextAlign(12);
        fit_legend->SetBorderSize(0);
        fit_legend->SetFillStyle(0);
        fit_legend->SetTextFont(42);
        fit_legend->SetTextSize(0.035);

        double mean_cs = p1 + 1.0;
        fit_legend->AddEntry((TObject*)nullptr, Form("Mean CS = %.3f #pm %.3f", mean_cs, ep1), "");

        fit_legend->Draw();

        std::string plot_title = obj ? obj->GetTitle() : "";
        drawATLASHeaderBlock(
            ndc_x0 - 0.03, ndc_y0 - 0.09,
            "Work in Progress",
            plot_title,
            32,
            kWhite, 0.70,
            kBlack, 0,
            0.01
        );
    }

    void styleRoI(TObject* obj, TCanvas* canvas, TClass* cl) {

    }

    void styleCSDistributionCombined(TObject* obj, TCanvas* canvas, TClass* cl) {

        auto stack = dynamic_cast<THStack*>(obj);
        if (!stack) return;

        stack->Draw("nostack hist");

        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(obj->GetTitle(), stack);

        const std::vector<Color_t> palette = {
            kAzure + 2,
            kGreen + 2,
            kOrange + 10,
            kMagenta + 2,
            kYellow - 3,
            kCyan - 4
        };
        TIter next(stack->GetHists());
        TH1* hist = nullptr;
        int index = 0;

        while ((hist = static_cast<TH1*>(next()))) {
            Color_t base_color = palette[index % palette.size()];

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

        drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.09,
            "Work in Progress",
            title_lines,
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

        TH1F* frame = canvas->DrawFrame(-9.0, 0.0, x_max, y_max);
        frame->SetTitle(h1->GetTitle());
        frame->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
        frame->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());

        applyATLASStyle(frame, canvas);

        h1->SetFillColorAlpha(kOrange - 2, 0.75);
        h1->SetFillStyle(1001);
        h1->SetLineColor(kOrange + 7);
        h1->SetLineWidth(2);
        h1->SetLineStyle(1);

        h1->Draw("HIST SAME");

        for (int i = 1; i < h1->GetNbinsX(); ++i) {
            double x = h1->GetXaxis()->GetBinUpEdge(i);

            double y_left = h1->GetBinContent(i);
            double y_right = h1->GetBinContent(i + 1);

            double y_max_line = std::min(y_left, y_right);

            if (y_max_line > 0) {
                TLine* edge = new TLine(x, 0.0, x, y_max_line);
                edge->SetLineColorAlpha(kOrange + 7, 0.55);
                edge->SetLineWidth(2);
                edge->Draw();
            }
        }

        TF1* fit = h1->GetFunction("gaus");
        if (!fit) fit = h1->GetFunction("half_gaus");

        if (!fit) {
            std::cerr << "No Gaussian fit found for histogram: " << h1->GetName() << std::endl;
            return;
        }

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

        canvas->RedrawAxis();

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TLegend* fit_legend = new TLegend(ndc_x0 + 0.03, ndc_y0 - 0.30, ndc_x0 + 0.35, ndc_y0 - 0.15);
        fit_legend->SetTextAlign(12);
        fit_legend->SetBorderSize(0);
        fit_legend->SetFillStyle(0);
        fit_legend->SetTextFont(42);
        fit_legend->SetTextSize(0.035);

        fit_legend->AddEntry(fit_band, "Gaussian Fit #pm 3#sigma Conf.", "fl");
        fit_legend->AddEntry((TObject*)nullptr, Form("#mu = %.3f #pm %.3f", p1, ep1), "");
        fit_legend->AddEntry((TObject*)nullptr, Form("#sigma = %.3f #pm %.3f", p2, ep2), "");

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

    void styleToFHeatmap(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto h2 = dynamic_cast<TH2*>(obj);
        if (!h2) return;

        double max_x = h2->GetXaxis()->GetXmax();
        h2->GetXaxis()->SetRangeUser(5200.0, max_x);

        h2->Draw("COLZ");

        auto [title_lines, x_label, y_label, legend_entries] = compilePlotLabels(h2->GetTitle(), h2);

        applyATLASStyle(obj, canvas);

        double ndc_x0 = canvas->GetLeftMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        drawATLASHeaderBlock(
            ndc_x0 + 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            title_lines,
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

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        std::string plot_title = obj ? obj->GetTitle() : "";
        drawATLASHeaderBlock(
            ndc_x0 - 0.03, ndc_y0 - 0.10,
            "Work in Progress",
            plot_title,
            32,
            kWhite, 0.00,
            kBlack, 0,
            0.01
        );
    }

    void styleToTCombSidesDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);

        auto [title_lines, x_label, y_label, compiled_legend_entries] = compilePlotLabels(obj->GetTitle(), stack);

        if (!stack) return;

        TList* hist_list = stack->GetHists();
        if (!hist_list || hist_list->GetSize() < 2) return;

        canvas->SetWindowSize(800, 840);
        canvas->SetCanvasSize(800, 840);

        canvas->cd();
        
        TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
        pad1->SetLeftMargin(0.16);
        pad1->SetRightMargin(0.05);
        pad1->SetTopMargin(0.07);
        pad1->SetBottomMargin(0.03);
        pad1->Draw();

        TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.30);
        pad2->SetLeftMargin(0.16);
        pad2->SetRightMargin(0.05);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.35);
        pad2->SetTicks(1, 1);
        pad2->Draw();

        std::vector<std::pair<double, Color_t>> mean_line_data;
        TH1* h1 = nullptr;
        TH1* h2 = nullptr;

        TIter next(hist_list);
        TH1* hist = nullptr;
        int index = 0;

        std::vector<TH1*> histograms;

        double global_percentile = 0.0;
        double percentile = 0.99;

        while ((hist = static_cast<TH1*>(next()))) {
            Color_t base_color = (index == 0) ? kAzure + 2 : kOrange + 10;
            mean_line_data.push_back({hist->GetMean(), base_color});

            if (index == 0) h1 = hist;
            if (index == 1) h2 = hist;

            hist->SetLineColor(base_color);
            hist->SetLineWidth(2);
            hist->SetLineStyle(1);

            Int_t trans_color = TColor::GetColorTransparent(base_color, 0.30);
            hist->SetFillColor(trans_color);
            hist->SetFillStyle(1001);

            if (hist->GetEffectiveEntries() > 0) {
                double q[1];
                double prob[1] = {percentile};
                hist->GetQuantiles(1, q, prob);
                if (q[0] > global_percentile) {
                    global_percentile = q[0];
                }
            }

            index++;
            histograms.push_back(hist);
        }

        pad1->cd();
        applyATLASStyle(obj, pad1);
        stack->Draw("nostack hist");

        if (stack->GetXaxis()) {
            stack->GetXaxis()->SetRangeUser(0.0, global_percentile);

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

        // GetUxmax() automatically reflects the new dynamic limit
        double x_max = pad1->GetUxmax();
        double y_min = pad1->GetUymin();
        double y_max = pad1->GetUymax();

        if (x_max > TOT_ROI_MAX + 1.0) {
            Objects::line(pad1, TOT_ROI_MAX, TOT_ROI_MAX, y_min, y_max, kBlack, 9, 1);
            Objects::hatchedRegion(pad1, TOT_ROI_MAX, x_max, y_min, y_max, 3244);

            double x_ndc = (TOT_ROI_MAX - pad1->GetX1()) / (pad1->GetX2() - pad1->GetX1());
            x_ndc += 0.015;
            double y_ndc = 0.35;

            TLatex* roi_text = new TLatex();
            roi_text->SetNDC(true);
            roi_text->SetTextAngle(90);
            roi_text->SetTextAlign(23);
            roi_text->SetTextFont(42);
            roi_text->SetTextSize(0.035);
            roi_text->SetTextColor(kBlack);
            roi_text->DrawLatex(x_ndc, y_ndc, "Streamer/Afterpulse Region");
        }

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
            ndc_y0 - 0.09,
            "Work in Progress",
            plot_title,
            32,
            kWhite, 0.70f,
            kWhite, 0,
            0.01
        );

        pad1->Modified();
        pad1->Update();

        std::vector<std::string> legend_entries = {" Side #eta_{1}", " Side #eta_{2}"};
        double legend_y = header ? header->GetY1NDC() - 0.05 : 0.70;
        double leg_width = 0.32;
        double leg_height = histograms.size() * 0.04;

        TLegend* leg = new TLegend(ndc_x0 - 0.05 - leg_width, legend_y - leg_height, ndc_x0 - 0.05, legend_y);
        leg->SetBorderSize(0);
        leg->SetLineColor(kWhite);
        leg->SetFillStyle(1001);
        leg->SetFillColorAlpha(kWhite, 0.5);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);

        for (size_t i = 0; i < histograms.size() && i < legend_entries.size(); ++i) {
            TH1* h = histograms[i];
            double mean_val = h->GetMean();
            std::string label_with_mean = Form("%s (#mu = %.1f ns)", legend_entries[i].c_str(), mean_val);

            leg->AddEntry(h, label_with_mean.c_str(), "f");
        }
        leg->Draw();

        pad2->cd();

        TH1* h_ratio = static_cast<TH1*>(h1->Clone("h_ratio"));
        h_ratio->SetDirectory(nullptr);
        h_ratio->Reset();
        h_ratio->Divide(h1, h2, 1.0, 1.0, "B");

        h_ratio->SetLineColor(kBlack);
        h_ratio->SetMarkerColor(kBlack);
        h_ratio->SetLineWidth(2);
        h_ratio->SetMarkerStyle(20);
        h_ratio->SetMarkerSize(0.8);
        h_ratio->SetFillStyle(0);

        h_ratio->Draw("ep"); 

        TAxis* rx = h_ratio->GetXaxis();
        TAxis* ry = h_ratio->GetYaxis();

        if (rx) {
            rx->SetRangeUser(0.0, global_percentile);
            rx->SetTitle("ToT [ns]");
            rx->SetTitleFont(42);
            rx->SetTitleSize(0.14);  
            rx->SetTitleOffset(1.0);
            rx->SetLabelFont(42);
            rx->SetLabelSize(0.11);
            rx->SetNdivisions(505);
            rx->SetTickLength(0.06);
        }

        if (ry) {
            ry->SetTitle("#eta_{1} / #eta_{2}");
            ry->SetTitleFont(42);
            ry->SetLabelFont(42);
            ry->SetTitleSize(0.14);
            ry->SetLabelSize(0.11);
            ry->SetTitleOffset(0.50);
            ry->SetNdivisions(505);
            h_ratio->SetMinimum(0.0); 
            h_ratio->SetMaximum(2.0);
        }

        pad2->Modified();
        pad2->Update();

        double x_max_pad2 = pad2->GetUxmax();
        double y_min_pad2 = pad2->GetUymin();
        double y_max_pad2 = pad2->GetUymax();

        Objects::line(pad2, 0.0, x_max_pad2, 1.0, 1.0, kGray+2, 2, 1);
        if (x_max_pad2 > TOT_ROI_MAX + 1.0) {
            Objects::line(pad2, TOT_ROI_MAX, TOT_ROI_MAX, y_min_pad2, y_max_pad2, kBlack, 9, 1);
            Objects::hatchedRegion(pad2, TOT_ROI_MAX, x_max_pad2, y_min_pad2, y_max_pad2, 3244);
        }

        pad2->Modified();
        pad2->Update();

        canvas->cd();
        canvas->Modified();
        canvas->Update();
    }

    void styleToTCombLayersDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);

        auto [title_lines, x_label, y_label, compiled_legend_entries] = compilePlotLabels(obj->GetTitle(), stack);

        if (!stack) return;

        TList* hist_list = stack->GetHists();
        if (!hist_list || hist_list->GetSize() < 3) return;

        canvas->SetWindowSize(800, 840);
        canvas->SetCanvasSize(800, 840);
        canvas->cd();

        TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
        pad1->SetLeftMargin(0.16);
        pad1->SetRightMargin(0.05);
        pad1->SetTopMargin(0.07);
        pad1->SetBottomMargin(0.03);
        pad1->Draw();

        TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.30);
        pad2->SetLeftMargin(0.16);
        pad2->SetRightMargin(0.05);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.35);
        pad2->SetTicks(1, 1);
        pad2->Draw();

        const std::vector<Color_t> palette = { kAzure + 2, kGreen + 2, kOrange + 10 };

        std::vector<std::pair<double, Color_t>> mean_line_data;
        std::vector<TH1*> histograms;

        TIter next(hist_list);
        TH1* hist = nullptr;
        int index = 0;

        double global_percentile = 0.0;
        double percentile = 0.99;

        while ((hist = static_cast<TH1*>(next()))) {
            Color_t base_color = palette[index];

            mean_line_data.push_back({hist->GetMean(), base_color});

            hist->SetLineColor(base_color);
            hist->SetLineWidth(2);
            hist->SetLineStyle(1);
            hist->SetFillStyle(0);

            // Calculate percentile for this histogram
            if (hist->GetEffectiveEntries() > 0) {
                double q[1];
                double prob[1] = {percentile};
                hist->GetQuantiles(1, q, prob);
                if (q[0] > global_percentile) {
                    global_percentile = q[0];
                }
            }

            histograms.push_back(hist);
            index++;
        }

        pad1->cd();
        applyATLASStyle(obj, pad1);
        stack->Draw("nostack hist");

        if (stack->GetXaxis()) {
            stack->GetXaxis()->SetRangeUser(0.0, global_percentile);

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

        if (x_max > TOT_ROI_MAX + 1.0) {
            Objects::line(pad1, TOT_ROI_MAX, TOT_ROI_MAX, y_min, y_max, kBlack, 9, 1);
            Objects::hatchedRegion(pad1, TOT_ROI_MAX, x_max, y_min, y_max, 3244);

            double x_ndc = (TOT_ROI_MAX - pad1->GetX1()) / (pad1->GetX2() - pad1->GetX1());
            x_ndc += 0.015;
            double y_ndc = 0.35;

            TLatex* roi_text = new TLatex();
            roi_text->SetNDC(true);
            roi_text->SetTextAngle(90);
            roi_text->SetTextAlign(23);
            roi_text->SetTextFont(42);
            roi_text->SetTextSize(0.035);
            roi_text->SetTextColor(kBlack);
            roi_text->DrawLatex(x_ndc, y_ndc, "Streamer/Afterpulse Region");
        }

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
            ndc_y0 - 0.09,
            "Work in Progress",
            plot_title,
            32,
            kWhite, 0.70f,
            kWhite, 0,
            0.01
        );

        pad1->Modified();
        pad1->Update();

        std::vector<std::string> legend_entries = {" Layer 0", " Layer 1", " Layer 2"};
        double legend_y = header ? header->GetY1NDC() - 0.10 : 0.70;
        double leg_width = 0.32;
        double leg_height = histograms.size() * 0.04;

        TLegend* leg = new TLegend(ndc_x0 - 0.05 - leg_width, legend_y - leg_height, ndc_x0 - 0.05, legend_y);
        leg->SetBorderSize(0);
        leg->SetLineColor(kWhite);
        leg->SetFillStyle(1001);
        leg->SetFillColorAlpha(kWhite, 0.5);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);

        for (size_t i = 0; i < histograms.size() && i < legend_entries.size(); ++i) {
            TH1* h = histograms[i];
            double mean_val = h->GetMean();
            std::string label_with_mean = Form("%s (#mu = %.1f ns)", legend_entries[i].c_str(), mean_val);
            leg->AddEntry(h, label_with_mean.c_str(), "l");
        }
        leg->Draw();

        pad2->cd();

        TH1* h_sum = static_cast<TH1*>(histograms[0]->Clone("h_sum"));
        h_sum->SetDirectory(nullptr);
        h_sum->Add(histograms[1]);
        h_sum->Add(histograms[2]);

        std::vector<TH1*> ratios;
        for (size_t i = 0; i < histograms.size(); ++i) {
            Color_t color = palette[i];

            TH1* h_ratio = static_cast<TH1*>(histograms[i]->Clone(Form("h_ratio_%zu", i)));
            h_ratio->SetDirectory(nullptr);
            h_ratio->Divide(histograms[i], h_sum, 1.0, 1.0, "B");

            h_ratio->SetLineColor(color);
            h_ratio->SetMarkerColor(color);
            h_ratio->SetLineWidth(2);
            h_ratio->SetMarkerStyle(20);
            h_ratio->SetMarkerSize(0.8);
            h_ratio->SetFillStyle(0);

            ratios.push_back(h_ratio);
        }

        ratios[0]->Draw("ep");

        TAxis* rx = ratios[0]->GetXaxis();
        TAxis* ry = ratios[0]->GetYaxis();

        if (rx) {
            rx->SetRangeUser(0.0, global_percentile);
            rx->SetTitle("ToT [ns]");
            rx->SetTitleFont(42);
            rx->SetTitleSize(0.14);
            rx->SetTitleOffset(1.0);
            rx->SetLabelFont(42);
            rx->SetLabelSize(0.11);
            rx->SetNdivisions(505);
            rx->SetTickLength(0.06);
        }

        if (ry) {
            ry->SetTitle("Layer / Total");
            ry->SetTitleFont(42);
            ry->SetLabelFont(42);
            ry->SetTitleSize(0.14);
            ry->SetLabelSize(0.11);
            ry->SetTitleOffset(0.50);
            ry->SetNdivisions(505);

            ratios[0]->SetMinimum(0.0);
            ratios[0]->SetMaximum(0.75);
        }

        ratios[1]->Draw("ep same");
        ratios[2]->Draw("ep same");

        pad2->Modified();
        pad2->Update();

        double x_max_pad2 = pad2->GetUxmax();
        double y_min_pad2 = pad2->GetUymin();
        double y_max_pad2 = pad2->GetUymax();

        Objects::line(pad2, 0.0, x_max_pad2, 1.0/3.0, 1.0/3.0, kGray+2, 2, 1);

        if (x_max_pad2 > TOT_ROI_MAX + 1.0) {
            Objects::line(pad2, TOT_ROI_MAX, TOT_ROI_MAX, y_min_pad2, y_max_pad2, kBlack, 9, 1);
            Objects::hatchedRegion(pad2, TOT_ROI_MAX, x_max_pad2, y_min_pad2, y_max_pad2, 3244);
        }

        pad2->Modified();
        pad2->Update();

        canvas->cd();
        canvas->Modified();
        canvas->Update();
    }

    void styleDelayDistribution(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);
        if (!stack) return;

        auto [title_lines, x_label, y_label, ignored_legend] = compilePlotLabels(obj->GetName(), stack);

        int idx = 0;
        TIter next(stack->GetHists());
        std::vector<TH1*> histograms;
        std::vector<std::string> custom_legend;
        TH1* hist = nullptr;

        double percentile = 0.99;
        double global_percentile = 0.0;

        std::vector<Color_t> pastel_colors = {kAzure + 1, kOrange - 2, kRed + 2};

        while ((hist = static_cast<TH1*>(next()))) {
            Color_t color = pastel_colors[idx % pastel_colors.size()];

            hist->SetLineColor(kBlack);
            hist->SetLineWidth(2);
            hist->SetLineStyle(1);

            hist->SetFillColor(color);
            hist->SetFillStyle(1001);

            if (hist->GetEffectiveEntries() > 0) {
                double q[1];
                double prob[1] = {percentile};
                hist->GetQuantiles(1, q, prob);
                if (q[0] > global_percentile) {
                    global_percentile = q[0];
                }
            }

            custom_legend.push_back(hist->GetTitle());
            idx++;
            histograms.push_back(hist);
        }

        stack->Draw("HIST");

        if (stack->GetXaxis()) {
            double absolute_max = stack->GetXaxis()->GetXmax();
            double dynamic_xMax = std::min(global_percentile + 2.0, absolute_max);

            stack->GetXaxis()->SetRangeUser(0.0, dynamic_xMax);
            stack->GetXaxis()->SetTitle(x_label.c_str());
        }
        if (stack->GetYaxis()) {
            stack->GetYaxis()->SetTitle(y_label.c_str());
        }

        applyATLASStyle(obj, canvas);
        enforceIntegerMinorTicks(stack->GetXaxis());
        canvas->Modified();
        canvas->Update();

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.03, ndc_y0 - 0.09,
            "Work in Progress", title_lines, 32,
            kWhite, 0.00,
            kBlack, 0, 0.01
        );

        canvas->Modified();
        canvas->Update();

        // Draw custom legend
        double legend_y = header ? header->GetY1NDC() - 0.02 : 0.70;
        double leg_width = 0.15;
        double leg_height = custom_legend.size() * 0.04;

        TLegend* leg = new TLegend(ndc_x0 - 0.03 - leg_width, legend_y - leg_height, ndc_x0 - 0.03, legend_y);
        leg->SetBorderSize(0);
        leg->SetLineColor(kWhite);
        leg->SetFillStyle(1001);
        leg->SetFillColorAlpha(kWhite, 0.5);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);

        for (size_t i = 0; i < histograms.size() && i < custom_legend.size(); ++i) {
            leg->AddEntry(histograms[i], custom_legend[i].c_str(), "f");
        }
        leg->Draw();

        canvas->RedrawAxis();
        canvas->Modified();
        canvas->Update();
    }

    void styleAvgDelayStrip(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);
        if (!stack) return;

        // Extract the title
        std::string full_title = obj->GetTitle();
        size_t semicolon_pos = full_title.find(';');
        std::string title = (semicolon_pos != std::string::npos) ? full_title.substr(0, semicolon_pos) : full_title;

        std::vector<TH1*> histograms;
        double global_max_with_err = -std::numeric_limits<double>::max();
        const std::vector<Color_t> palette = {
            kAzure + 2,
            kGreen + 2,
            kOrange + 10,
            kMagenta + 2,
            kYellow - 3,
            kCyan - 4
        };
        if (TList* hists = stack->GetHists()) {
            TIter next(hists);
            TObject* hist_obj;
            int idx = 0;

            while ((hist_obj = next())) {
                if (auto h = dynamic_cast<TH1*>(hist_obj)) {
                    Color_t color = palette[idx % palette.size()];

                    h->SetMarkerStyle(8);
                    h->SetMarkerSize(1);
                    h->SetMarkerColor(color);
                    h->SetLineColor(color);
                    h->SetLineWidth(1);

                    // Extract max + error for this specific histogram
                    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
                        double val_with_err = h->GetBinContent(bin) + h->GetBinError(bin);
                        if (val_with_err > global_max_with_err) {
                            global_max_with_err = val_with_err;
                        }
                    }

                    histograms.push_back(h);
                    idx++;
                }
            }
        }

        stack->Draw("NOSTACK PE");

        if (auto h = stack->GetHistogram()) {
            if (TAxis* yAxis = h->GetYaxis()) {
                if (global_max_with_err <= 0) global_max_with_err = 1.0;
                yAxis->SetRangeUser(0.0, global_max_with_err * 1.25);
            }
        }

        if (auto named_obj = dynamic_cast<TNamed*>(obj)) named_obj->SetTitle(title.c_str());

        // Update canvas offsets, margins and enforce integer minor ticks
        applyATLASStyle(obj, canvas);
        if (stack->GetHistogram()) enforceIntegerMinorTicks(stack->GetHistogram()->GetXaxis());
        canvas->Modified();
        canvas->Update();

        // Draw the ATLAS header block
        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.03, ndc_y0 - 0.09,
            "Work in Progress", title, 32,
            kWhite, 0.70, kBlack, 0, 0.01
        );

        canvas->Modified();
        canvas->Update();

        // Draw the legend
        std::vector<std::string> legend_entries = {"Layer 0", "Layer 1", "Layer 2"};
        double legend_y = header ? header->GetY1NDC() - 0.02 : 0.70;
        double leg_width = 0.15;
        double leg_height = histograms.size() * 0.04;

        TLegend* leg = new TLegend(ndc_x0 - 0.03 - leg_width, legend_y - leg_height, ndc_x0 - 0.03, legend_y);
        leg->SetBorderSize(0);
        leg->SetLineColor(kWhite);
        leg->SetFillStyle(1001);
        leg->SetFillColorAlpha(kWhite, 0.5);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);

        for (size_t i = 0; i < histograms.size() && i < legend_entries.size(); ++i) {
            leg->AddEntry(histograms[i], legend_entries[i].c_str(), "pe");
        }
        leg->Draw();

        canvas->Modified();
        canvas->Update();
    }

    void styleAvgMultStrip(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);
        if (!stack) return;

        // Extract the title
        std::string full_title = obj->GetTitle();
        size_t semicolon_pos = full_title.find(';');
        std::string title = (semicolon_pos != std::string::npos) ? full_title.substr(0, semicolon_pos) : full_title;

        // Draw the stack with points and error bars
        stack->Draw("NOSTACK PE");

        if (auto h = stack->GetHistogram()) {
            if (TAxis* yAxis = h->GetYaxis()) {
                setRange(stack, yAxis, AxisType::Y, 0.95, std::nullopt, {.y_max = 2.0});
            }
        }

        std::vector<TH1*> histograms;
        const std::vector<Color_t> palette = {
            kAzure + 2,
            kGreen + 2,
            kOrange + 10,
            kMagenta + 2,
            kYellow - 3,
            kCyan - 4
        };
        if (TList* hists = stack->GetHists()) {
            TIter next(hists);
            TObject* hist_obj;
            int idx = 0;

            while ((hist_obj = next())) {
                if (auto h = dynamic_cast<TH1*>(hist_obj)) {
                    Color_t color = palette[idx % palette.size()];

                    h->SetMarkerStyle(8);
                    h->SetMarkerSize(1);
                    h->SetMarkerColor(color);
                    h->SetLineColor(color);
                    h->SetLineWidth(1);

                    idx++;
                    histograms.push_back(h);
                }
            }
        }

        // Draw a horizonal line at y=1.0 to indicate the baseline
        double x_min = stack->GetXaxis()->GetXmin();
        double x_max = stack->GetXaxis()->GetXmax();
        TLine* baseline_line = new TLine(x_min, 1.0, x_max, 1.0);
        baseline_line->SetLineColor(kBlack);
        baseline_line->SetLineStyle(2);
        baseline_line->SetLineWidth(1);
        baseline_line->Draw();

        if (auto named_obj = dynamic_cast<TNamed*>(obj)) named_obj->SetTitle(title.c_str());

        applyATLASStyle(obj, canvas);
        enforceIntegerMinorTicks(stack->GetXaxis());
        canvas->Modified();
        canvas->Update();

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.03, ndc_y0 - 0.09,
            "Work in Progress", title, 32,
            kWhite, 0.70, kBlack, 0, 0.01
        );

        canvas->Modified();
        canvas->Update();

        // Draw legend
        std::vector<std::string> legend_entries = {"Layer 0", "Layer 1", "Layer 2"};
        double legend_y = header ? header->GetY1NDC() - 0.02 : 0.70;
        double leg_width = 0.15;
        double leg_height = histograms.size() * 0.04;

        TLegend* leg = new TLegend(ndc_x0 - 0.03 - leg_width, legend_y - leg_height, ndc_x0 - 0.03, legend_y);
        leg->SetBorderSize(0);
        leg->SetLineColor(kWhite);
        leg->SetFillStyle(1001);
        leg->SetFillColorAlpha(kWhite, 0.5);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);

        for (size_t i = 0; i < histograms.size() && i < legend_entries.size(); ++i) {
            leg->AddEntry(histograms[i], legend_entries[i].c_str(), "pe");
        }
        leg->Draw();

        canvas->Modified();
        canvas->Update();
    }

    void styleFracMultStrip(TObject* obj, TCanvas* canvas, TClass* cl) {
        auto stack = dynamic_cast<THStack*>(obj);
        if (!stack) return;

        // Extract the title
        std::string full_title = obj->GetTitle();
        size_t semicolon_pos = full_title.find(';');
        std::string title = (semicolon_pos != std::string::npos) ? full_title.substr(0, semicolon_pos) : full_title;

        std::vector<TH1*> histograms;
        const std::vector<Color_t> palette = {
            kAzure + 2,
            kGreen + 2,
            kOrange + 10,
            kMagenta + 2,
            kYellow - 3,
            kCyan - 4
        };
        if (TList* hists = stack->GetHists()) {
            TIter next(hists);
            TObject* hist_obj;
            int idx = 0;

            while ((hist_obj = next())) {
                if (auto h = dynamic_cast<TH1*>(hist_obj)) {
                    Color_t color = palette[idx % palette.size()];

                    h->SetMarkerStyle(8);
                    h->SetMarkerSize(1);
                    h->SetMarkerColor(color);
                    h->SetLineColor(color);
                    h->SetLineWidth(1);

                    // Hide markers/errors for empty bins so they don't clutter the baseline
                    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
                        if (h->GetBinContent(bin) == 0) {
                            h->SetBinError(bin, 0.0);
                        }
                    }

                    idx++;
                    histograms.push_back(h);
                }
            }
        }

        stack->Draw("NOSTACK PE");

        applyATLASStyle(obj, canvas);
        enforceIntegerMinorTicks(stack->GetXaxis());
        canvas->Modified();
        canvas->Update();

        double ndc_x0 = 1.0 - canvas->GetRightMargin();
        double ndc_y0 = 1.0 - canvas->GetTopMargin();

        TPaveText* header = drawATLASHeaderBlock(
            ndc_x0 - 0.03, ndc_y0 - 0.09,
            "Work in Progress", title, 32,
            kWhite, 0.70, kBlack, 0, 0.01
        );

        canvas->Modified();
        canvas->Update();

        // Draw legend
        std::vector<std::string> legend_entries = {"Layer 0", "Layer 1", "Layer 2"};
        double legend_y = header ? header->GetY1NDC() - 0.02 : 0.70;
        double leg_width = 0.15;
        double leg_height = histograms.size() * 0.04;

        TLegend* leg = new TLegend(ndc_x0 - 0.03 - leg_width, legend_y - leg_height, ndc_x0 - 0.03, legend_y);
        leg->SetBorderSize(0);
        leg->SetLineColor(kWhite);
        leg->SetFillStyle(1001);
        leg->SetFillColorAlpha(kWhite, 0.5);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);

        for (size_t i = 0; i < histograms.size() && i < legend_entries.size(); ++i) {
            leg->AddEntry(histograms[i], legend_entries[i].c_str(), "pe");
        }
        leg->Draw();

        canvas->Modified();
        canvas->Update();
    }

    void styleDefaultPlot(TObject* obj, TCanvas* canvas, TClass* cl) {

        if (cl->InheritsFrom(TH2::Class())) {
            obj->Draw("COLZ");
        } else if (cl->InheritsFrom(TH1::Class())) {
            obj->Draw("HIST");
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