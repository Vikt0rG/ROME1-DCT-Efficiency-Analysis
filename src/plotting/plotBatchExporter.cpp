#include <iostream>
#include <regex>
#include <filesystem>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TLegend.h>

#include "plotting/plotStyler.hpp"
#include "plotting/plotBatchExporter.hpp"
#include "core/constants.hpp"

namespace PlotterHelpers {
namespace BatchExporter {

namespace {
    // Unnamed namespace keeps scanDirectory locally scoped to this compilation block
    void scanDirectory(TDirectory* dir, const std::filesystem::path& current_output_path) {
        TIter next(dir->GetListOfKeys());
        TKey* key = nullptr;

        while ((key = static_cast<TKey*>(next()))) {
            TClass* cl = TClass::GetClass(key->GetClassName());
            if (!cl) continue;

            if (cl->InheritsFrom(TDirectory::Class())) {
                TDirectory* subDir = dynamic_cast<TDirectory*>(key->ReadObj());
                if (subDir) {
                    std::string sub_dir_name = key->GetName();
                    std::filesystem::path next_output_path;
                    if (sub_dir_name == "analysis") {
                        next_output_path = current_output_path;
                    } else {
                        next_output_path = current_output_path / sub_dir_name;
                        std::filesystem::create_directories(next_output_path);
                    }
                    std::cout << "[ATLAS Export] Entering subdirectory: " << sub_dir_name << std::endl;
                    std::filesystem::create_directories(next_output_path);
                    scanDirectory(subDir, next_output_path);
                }
                continue;
            }

            if (cl->InheritsFrom(TH1::Class()) ||
                cl->InheritsFrom(TGraph::Class()) ||
                cl->InheritsFrom(TMultiGraph::Class()) ||
                cl->InheritsFrom(TGraphAsymmErrors::Class()) ||
                cl->InheritsFrom(THStack::Class())) {
                TObject* obj = key->ReadObj();
                if (!obj) continue;

                TCanvas* canvas = new TCanvas("c_temp", "", 800, 600);
                canvas->cd();

                std::string obj_name = obj->GetName();
                
                // Call across module boundary into our clean cosmetic classification tables
                PlotCategory category = PlotStyler::getPlotCategory(obj);
                auto custom_styler = PlotStyler::getCustomStyler(category);

                if (custom_styler) {
                    custom_styler(obj, canvas, cl);
                } else {
                    PlotStyler::styleDefaultPlot(obj, canvas, cl);
                }

                std::filesystem::path export_path = current_output_path / (obj_name + ".pdf");
                canvas->SaveAs(export_path.string().c_str());
                
                delete canvas;
                delete obj;
            }
        }
    }
} // namespace

void autoExportToATLASPDF(const std::string& root_file_path, const std::filesystem::path& target_plots_dir) {
    gSystem->mkdir(target_plots_dir.string().c_str(), kTRUE);
    TFile* file = TFile::Open(root_file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        if (file) delete file;
        throw std::runtime_error("Export pipeline failed to open input file: " + root_file_path);
    }

    std::cout << "[ATLAS Export] Beginning recursive sweep of " << root_file_path << "..." << std::endl;
    scanDirectory(file, target_plots_dir);
    std::cout << "[ATLAS Export] Completed cleanly. Plots saved to: " << target_plots_dir << std::endl;

    file->Close();
    delete file;
}

// Helper function to build a global TMultiGraph for global RPC plots across all layers for a given metric
void buildGlobalMultiGraphs(TDirectory* config_dir, const std::filesystem::path& config_output_path) {
    if (!config_dir) return;

    // Dynamically find any available group directory to act as a "blueprint"
    TDirectory* blueprint_dir = nullptr;
    TIter next_top_key(config_dir->GetListOfKeys());
    TKey* top_key = nullptr;

    while ((top_key = static_cast<TKey*>(next_top_key()))) {
        TClass* cl = TClass::GetClass(top_key->GetClassName());
        if (cl && cl->InheritsFrom(TDirectory::Class())) {
            std::string dir_name = top_key->GetName();
            if (dir_name.rfind("group_", 0) == 0) {
                blueprint_dir = dynamic_cast<TDirectory*>(top_key->ReadObj());
                break; 
            }
        }
    }

    if (!blueprint_dir) {
        std::cerr << "Error: No group_X directories found in the configuration directory. Cannot build global multi-graphs." << std::endl;
        return;
    }

    // Iterate through analysis categories (e.g., efficiency_analysis)
    TIter next_analysis_dir(blueprint_dir->GetListOfKeys());
    TKey* analysis_key = nullptr;

    while ((analysis_key = static_cast<TKey*>(next_analysis_dir()))) {
        TClass* cl_dir = TClass::GetClass(analysis_key->GetClassName());
        if (!cl_dir || !cl_dir->InheritsFrom(TDirectory::Class())) continue;

        std::string analysis_subdir_name = analysis_key->GetName();
        TDirectory* analysis_dir = dynamic_cast<TDirectory*>(analysis_key->ReadObj());
        if (!analysis_dir) continue;

        // Differentiate between Strip Metrics and Layer Metrics
        bool is_strip_analysis = (analysis_subdir_name == "tot_analysis" ||
                                  analysis_subdir_name == "multiplicity_analysis" ||
                                  analysis_subdir_name == "noise_rate_strips_analysis");

        if (is_strip_analysis) {
            // STRIP METRIC LOGIC: Fetch from layerX folders inside each group and assemble MGs
            TIter next_group(config_dir->GetListOfKeys());
            TKey* group_key = nullptr;

            while ((group_key = static_cast<TKey*>(next_group()))) {
                std::string group_dir_name = group_key->GetName();
                if (group_dir_name.rfind("group_", 0) != 0) continue;
                std::string prefix = "group_";
                std::string clean_group = group_dir_name.substr(prefix.length());

                for (int layer = 0; layer < LAYER_COUNT; ++layer) {
                    std::string layer_folder = "layer" + std::to_string(layer);
                    std::string target_path = group_dir_name + "/" + analysis_subdir_name + "/" + layer_folder;

                    TDirectory* layer_dir = config_dir->GetDirectory(target_path.c_str());
                    if (!layer_dir) continue;

                    std::map<std::string, TMultiGraph*> mg_map;
                    TIter next_strip(layer_dir->GetListOfKeys());
                    TKey* strip_key = nullptr;

                    while ((strip_key = static_cast<TKey*>(next_strip()))) {
                        TObject* obj = strip_key->ReadObj();
                        if (auto g = dynamic_cast<TGraph*>(obj)) {

                            // Check if this layer was actually scanned by looking at its X values
                            double x_start = 0.0, x_middle = 0.0, x_end = 0.0, y_dummy = 0.0;
                            g->GetPoint(0, x_start, y_dummy);
                            g->GetPoint(g->GetN() / 2, x_middle, y_dummy);
                            g->GetPoint(g->GetN() - 1, x_end, y_dummy);

                            // Skip graphs that have no meaningful scan range (e.g., flatlined layers)
                            if (std::abs(x_end - x_start) < 1.0 && std::abs(x_middle - x_start) < 1.0) {
                                delete obj;
                                break;
                            }

                            std::string g_name = g->GetName(); 

                            std::string base_metric = g_name;
                            size_t strip_pos = g_name.rfind("_strip");
                            if (strip_pos != std::string::npos) {
                                base_metric = g_name.substr(0, strip_pos);
                            }

                            // Avoid names like "layer0_layer0"
                            std::string dynamic_suffix = "_" + clean_group;
                            if (clean_group != layer_folder) {
                                dynamic_suffix += "_" + layer_folder;
                            }

                            if (mg_map.find(base_metric) == mg_map.end()) {
                                TMultiGraph* strip_mg = new TMultiGraph();
                                std::string mg_name = base_metric + dynamic_suffix;
                                strip_mg->SetName(mg_name.c_str());
                                strip_mg->SetTitle((mg_name + ";High Voltage [V];Value").c_str());
                                mg_map[base_metric] = strip_mg;
                            }

                            TGraph* clone = static_cast<TGraph*>(g->Clone());
                            clone->SetName(g->GetName());
                            clone->SetTitle(g->GetTitle());
                            mg_map[base_metric]->Add(clone, "P");
                        }
                        delete obj;
                    }

                    for (auto const& [base_metric, strip_mg] : mg_map) {

                        TCanvas* canvas = new TCanvas("c", "", 800, 600);
                        canvas->cd();

                        PlotCategory category = PlotterHelpers::PlotStyler::getPlotCategory(strip_mg);
                        auto custom_styler = PlotterHelpers::PlotStyler::getCustomStyler(category);

                        if (custom_styler) custom_styler(strip_mg, canvas, TMultiGraph::Class());
                        else PlotterHelpers::PlotStyler::styleDefaultPlot(strip_mg, canvas, TMultiGraph::Class());

                        // Use the same smart naming logic for the PDF output
                        std::string dynamic_suffix = "_" + clean_group;
                        if (clean_group != layer_folder) {
                            dynamic_suffix += "_" + layer_folder;
                        }

                        std::string final_mg_name = base_metric + dynamic_suffix;
                        std::filesystem::path export_file = config_output_path / analysis_subdir_name / (final_mg_name + ".pdf");
                        std::filesystem::create_directories(export_file.parent_path());

                        canvas->SaveAs(export_file.string().c_str());
                        delete canvas;
                        delete strip_mg;
                    }
                }   // End of layer loop
            }   // End of group loop
        } else {
            // LAYER & GLOBAL METRIC LOGIC: Handle 1D stitching and 2D Heatmaps
            TIter next_metric_key(analysis_dir->GetListOfKeys());
            TKey* metric_key = nullptr;

            while ((metric_key = static_cast<TKey*>(next_metric_key()))) {
                TClass* cl_metric = TClass::GetClass(metric_key->GetClassName());
                if (!cl_metric) continue; 

                std::string metric_name = metric_key->GetName();
                TObject* blueprint_obj = metric_key->ReadObj();
                if (blueprint_obj->InheritsFrom(TH1::Class())) {
                    static_cast<TH1*>(blueprint_obj)->SetDirectory(nullptr);
                }

                std::vector<std::pair<TObject*, std::string>> export_queue;
                // ------------------------------------------------------------------
                // CASE A: 2D Heatmaps
                if (cl_metric->InheritsFrom(TH2::Class())) {
                    TIter next_group(config_dir->GetListOfKeys());
                    TKey* group_key = nullptr;

                    while ((group_key = static_cast<TKey*>(next_group()))) {
                        std::string prefix = "group_";
                        std::string group_dir_name = group_key->GetName();
                        if (group_dir_name.rfind(prefix, 0) != 0) continue;

                        std::string clean_group = group_dir_name.substr(prefix.length());
                        std::string target_path = group_dir_name + "/" + analysis_subdir_name + "/" + metric_name;

                        TH2* raw_h2 = config_dir->Get<TH2>(target_path.c_str());
                        if (!raw_h2) continue;

                        TH2* h2_clone = static_cast<TH2*>(raw_h2->Clone());
                        h2_clone->SetDirectory(nullptr);
                        h2_clone->SetTitle((std::string(h2_clone->GetTitle()) + "_" + clean_group).c_str());

                        export_queue.push_back({h2_clone, metric_name + "_" + clean_group});
                    }
                }
                // ------------------------------------------------------------------
                // CASE B: 1D Graphs
                else if (cl_metric->InheritsFrom(TGraph::Class()) || cl_metric->InheritsFrom(TMultiGraph::Class())) {
                    TMultiGraph* global_mg = new TMultiGraph();
                    global_mg->SetName(metric_name.c_str());
                    if (auto named = dynamic_cast<TNamed*>(blueprint_obj)) {
                        global_mg->SetTitle(named->GetTitle());
                    }

                    bool graph_added = false;

                    TIter next_group(config_dir->GetListOfKeys());
                    TKey* group_key = nullptr;
                    while ((group_key = static_cast<TKey*>(next_group()))) {
                        std::string prefix = "group_";
                        std::string group_dir_name = group_key->GetName();
                        if (group_dir_name.rfind(prefix, 0) != 0) continue;

                        std::string clean_group = group_dir_name.substr(prefix.length());
                        std::string target_path = group_dir_name + "/" + analysis_subdir_name + "/" + metric_name;

                        TObject* scan_obj = config_dir->Get(target_path.c_str());
                        if (!scan_obj) continue;

                        if (auto mg = dynamic_cast<TMultiGraph*>(scan_obj)) {
                            if (mg->GetListOfGraphs()) {
                                for (TObject* gr_obj : *mg->GetListOfGraphs()) {
                                    auto g = static_cast<TGraph*>(gr_obj);

                                    double x_start = 0.0, x_middle = 0.0, x_end = 0.0, y_dummy = 0.0;
                                    g->GetPoint(0, x_start, y_dummy);
                                    g->GetPoint(g->GetN() / 2, x_middle, y_dummy);
                                    g->GetPoint(g->GetN() - 1, x_end, y_dummy);

                                    if (std::abs(x_end - x_start) < 1.0 && std::abs(x_middle - x_start) < 1.0) continue;

                                    std::string g_name = g->GetName();
                                    TGraph* clone = static_cast<TGraph*>(g->Clone());
                                    clone->SetName((metric_name + "_" + clean_group + "_" + g_name).c_str());

                                    if (clean_group.find("layer") == 0) clone->SetTitle(g->GetTitle());
                                    else clone->SetTitle((clean_group + " - " + g->GetTitle()).c_str()); 

                                    global_mg->Add(clone, "P");
                                    graph_added = true;
                                }
                            }
                        } else if (auto g = dynamic_cast<TGraph*>(scan_obj)) {
                            TGraph* clone = static_cast<TGraph*>(g->Clone());
                            clone->SetName((metric_name + "_" + clean_group).c_str());
                            clone->SetTitle(clean_group.c_str());
                            global_mg->Add(clone, "P");
                            graph_added = true;
                        }
                        delete scan_obj; 
                    }
                    
                    if (graph_added) {
                        export_queue.push_back({global_mg, metric_name});
                    } else {
                        delete global_mg;
                    }
                }

                // ------------------------------------------------------------------
                // Export all queued objects to PDF
                for (auto& [obj_to_export, filename_base] : export_queue) {
                    TCanvas* canvas = new TCanvas("c", "", 800, 600);
                    canvas->cd();

                    TClass* obj_class = obj_to_export->IsA();

                    PlotCategory category = PlotterHelpers::PlotStyler::getPlotCategory(obj_to_export);
                    auto custom_styler = PlotterHelpers::PlotStyler::getCustomStyler(category);

                    if (custom_styler) custom_styler(obj_to_export, canvas, obj_class);
                    else PlotterHelpers::PlotStyler::styleDefaultPlot(obj_to_export, canvas, obj_class);

                    std::filesystem::path export_file = config_output_path / analysis_subdir_name / (filename_base + ".pdf");
                    std::filesystem::create_directories(export_file.parent_path());

                    canvas->SaveAs(export_file.string().c_str());

                    delete canvas;
                    delete obj_to_export; 
                }

                delete blueprint_obj;
            }
        }
        delete analysis_dir;
    }
    delete blueprint_dir;

    TDirectory* summary_dir = config_dir->GetDirectory("cross_group_summaries");
    if (!summary_dir) return;

    TIter next_summary(summary_dir->GetListOfKeys());
    TKey* summary_key = nullptr;

    while ((summary_key = static_cast<TKey*>(next_summary()))) {
        TObject* obj = summary_key->ReadObj();
        if (!obj) continue;

        if (!obj->InheritsFrom(TGraph::Class())) {
            delete obj;
            continue;
        }

        TCanvas* canvas = new TCanvas("c_summary", "", 800, 600);
        canvas->cd();

        PlotCategory category = PlotterHelpers::PlotStyler::getPlotCategory(obj);
        auto custom_styler = PlotterHelpers::PlotStyler::getCustomStyler(category);

        if (custom_styler) custom_styler(obj, canvas, obj->IsA());
        else PlotterHelpers::PlotStyler::styleDefaultPlot(obj, canvas, obj->IsA());

        // Save into a dedicated 'cross_group_summaries' folder in the PDF output directory
        std::filesystem::path export_file = config_output_path / "cross_group_summaries" / (std::string(obj->GetName()) + ".pdf");
        std::filesystem::create_directories(export_file.parent_path());

        canvas->SaveAs(export_file.string().c_str());

        delete canvas;
        delete obj;
    }
}

} // namespace BatchExporter
} // namespace PlotterHelpers