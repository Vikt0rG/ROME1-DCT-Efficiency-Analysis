#include <iostream>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TH1.h>
#include <THStack.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TSystem.h>

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

        // STEP 1: Dynamically find any available scanned_layer directory to act as a "blueprint"
        TDirectory* blueprint_dir = nullptr;
        TIter next_top_key(config_dir->GetListOfKeys());
        TKey* top_key = nullptr;

        while ((top_key = static_cast<TKey*>(next_top_key()))) {
            TClass* cl = TClass::GetClass(top_key->GetClassName());
            if (cl && cl->InheritsFrom(TDirectory::Class())) {
                std::string dir_name = top_key->GetName();
                if (dir_name.rfind("scanned_layer_", 0) == 0) {
                    blueprint_dir = dynamic_cast<TDirectory*>(top_key->ReadObj());
                    break; 
                }
            }
        }

        if (!blueprint_dir) {
            std::cerr << "Error: No scanned_layer_X directories found in the configuration directory. Cannot build global multi-graphs." << std::endl;
            return;
        }

        // STEP 2: Iterate through analysis categories (e.g., efficiency_analysis)
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
                // STRIP METRIC LOGIC: Fetch directly from layerX folders and assemble MGs per metric
                for (int layer = 0; layer < LAYER_COUNT; ++layer) {
                    std::string scan_folder = "scanned_layer_" + std::to_string(layer);
                    std::string layer_folder = "layer" + std::to_string(layer);
                    std::string target_path = scan_folder + "/" + analysis_subdir_name + "/" + layer_folder;

                    TDirectory* layer_dir = config_dir->GetDirectory(target_path.c_str());
                    if (!layer_dir) continue; // Safely skip if this scan or layer folder is missing

                    // Map to hold distinct MultiGraphs for different metrics (eta1, eta2, etc.)
                    std::map<std::string, TMultiGraph*> mg_map;

                    TIter next_strip(layer_dir->GetListOfKeys());
                    TKey* strip_key = nullptr;

                    // Grab all single strip graphs inside the layerX folder and sort them
                    while ((strip_key = static_cast<TKey*>(next_strip()))) {
                        TObject* obj = strip_key->ReadObj();
                        if (auto g = dynamic_cast<TGraph*>(obj)) {
                            std::string g_name = g->GetName(); // e.g., "avg_tot_eta1_strip0"
                            
                            // Extract the base metric name (e.g., "avg_tot_eta1")
                            std::string base_metric = g_name;
                            size_t strip_pos = g_name.rfind("_strip");
                            if (strip_pos != std::string::npos) {
                                base_metric = g_name.substr(0, strip_pos);
                            }

                            // Initialize a new TMultiGraph for this specific metric if it doesn't exist yet
                            if (mg_map.find(base_metric) == mg_map.end()) {
                                TMultiGraph* strip_mg = new TMultiGraph();
                                std::string mg_name = base_metric + "_layer" + std::to_string(layer);
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

                    // Now style and save all the distinct MultiGraphs we just built
                    for (auto const& [base_metric, strip_mg] : mg_map) {
                        TCanvas* canvas = new TCanvas("c", "", 800, 600);
                        canvas->cd();

                        PlotCategory category = PlotterHelpers::PlotStyler::getPlotCategory(strip_mg);
                        auto custom_styler = PlotterHelpers::PlotStyler::getCustomStyler(category);

                        if (custom_styler) custom_styler(strip_mg, canvas, TMultiGraph::Class());
                        else PlotterHelpers::PlotStyler::styleDefaultPlot(strip_mg, canvas, TMultiGraph::Class());

                        // Save nicely as e.g., config_output/tot_analysis/avg_tot_eta1_layer0.pdf
                        std::string final_mg_name = base_metric + "_layer" + std::to_string(layer);
                        std::filesystem::path export_file = config_output_path / analysis_subdir_name / (final_mg_name + ".pdf");
                        std::filesystem::create_directories(export_file.parent_path());
                        
                        canvas->SaveAs(export_file.string().c_str());

                        delete canvas;
                        delete strip_mg;
                    }
                }
            }
            else {
                // LAYER METRIC LOGIC: Stitch layer0, layer1, layer2 across scans
                TIter next_metric_key(analysis_dir->GetListOfKeys());
                TKey* metric_key = nullptr;

                while ((metric_key = static_cast<TKey*>(next_metric_key()))) {
                    TClass* cl_metric = TClass::GetClass(metric_key->GetClassName());

                    // IGNORE the layer0/layer1 folders. Only process top-level MGs/Graphs.
                    if (!cl_metric || !(cl_metric->InheritsFrom(TGraph::Class()) || cl_metric->InheritsFrom(TMultiGraph::Class()))) {
                        continue; 
                    }

                    std::string metric_name = metric_key->GetName();
                    TObject* blueprint_obj = metric_key->ReadObj();

                    // Create the new Global MultiGraph to hold our 3 stitched layers
                    TMultiGraph* global_mg = new TMultiGraph();
                    global_mg->SetName(metric_name.c_str());

                    if (auto named = dynamic_cast<TNamed*>(blueprint_obj)) {
                        global_mg->SetTitle(named->GetTitle());
                    }

                    bool graph_added = false;

                    // Fetch Layer 0 from Scan 0, Layer 1 from Scan 1, Layer 2 from Scan 2
                    for (int scan = 0; scan < 3; ++scan) {
                        std::string target_path = "scanned_layer_" + std::to_string(scan) + "/" + 
                                                analysis_subdir_name + "/" + metric_name;

                        TObject* scan_obj = config_dir->Get(target_path.c_str());
                        if (!scan_obj) continue;

                        TGraph* extracted_graph = nullptr;

                        if (auto mg = dynamic_cast<TMultiGraph*>(scan_obj)) {
                            if (mg->GetListOfGraphs()) {
                                int idx = 0;
                                for (TObject* gr_obj : *mg->GetListOfGraphs()) {
                                    auto g = static_cast<TGraph*>(gr_obj);
                                    std::string g_name = g->GetName();

                                    // Find the graph corresponding to the current scan layer
                                    if (g_name.find("layer" + std::to_string(scan)) != std::string::npos) {
                                        extracted_graph = g;
                                        break;
                                    } else if (idx == scan) {
                                        extracted_graph = g; // Fallback to index
                                    }
                                    idx++;
                                }
                            }
                        } else if (auto g = dynamic_cast<TGraph*>(scan_obj)) {
                            extracted_graph = g; // If it's already a flat graph
                        }

                        // Clone it and add it to the global multigraph
                        if (extracted_graph) {
                            TGraph* clone = static_cast<TGraph*>(extracted_graph->Clone());
                            clone->SetName((metric_name + "_layer" + std::to_string(scan)).c_str());
                            global_mg->Add(clone, "P");
                            graph_added = true;
                        }
                        
                        delete scan_obj; 
                    }

                    // Style and Render the Global Plot
                    if (graph_added) {
                        TCanvas* canvas = new TCanvas("c", "", 800, 600);
                        canvas->cd();

                        PlotCategory category = PlotterHelpers::PlotStyler::getPlotCategory(global_mg);
                        auto custom_styler = PlotterHelpers::PlotStyler::getCustomStyler(category);

                        if (custom_styler) {
                            custom_styler(global_mg, canvas, TMultiGraph::Class());
                        } else {
                            PlotterHelpers::PlotStyler::styleDefaultPlot(global_mg, canvas, TMultiGraph::Class());
                        }

                        // Export to config_output/analysis_folder/metric_name.pdf
                        std::filesystem::path export_file = config_output_path / analysis_subdir_name / (metric_name + ".pdf");
                        std::filesystem::create_directories(export_file.parent_path());

                        canvas->SaveAs(export_file.string().c_str());
                        delete canvas;
                    }

                    delete global_mg;
                    delete blueprint_obj;
                }
            }
            delete analysis_dir;
        }
        delete blueprint_dir;
    }

} // namespace BatchExporter
} // namespace PlotterHelpers