#include <iostream>

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TSystem.h>

#include "plotting/plotStyler.hpp"
#include "plotting/plotBatchExporter.hpp"

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
                    std::filesystem::path next_output_path = current_output_path / sub_dir_name;
                    std::cout << "[ATLAS Export] Entering subdirectory: " << sub_dir_name << std::endl;
                    std::filesystem::create_directories(next_output_path);
                    scanDirectory(subDir, next_output_path);
                }
                continue;
            }

            if (cl->InheritsFrom(TH1::Class()) || cl->InheritsFrom(TGraph::Class()) || cl->InheritsFrom(TMultiGraph::Class())) {
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

        // STEP 1: Dynamically find any available scanned_layer directory to act as a blueprint
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
            std::cerr << "Warning: No 'scanned_layer_*' directories found inside configuration: " 
                    << config_dir->GetName() << ". Skipping." << std::endl;
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

            // STEP 3: Dynamically open "layer0" (or layer1/layer2 if 0 is missing) to extract real metric names
            TDirectory* layer_blueprint = nullptr;
            for (int l = 0; l < 3; ++l) {
                layer_blueprint = dynamic_cast<TDirectory*>(analysis_dir->Get(("layer" + std::to_string(l)).c_str()));
                if (layer_blueprint) break;
            }

            if (!layer_blueprint) {
                delete analysis_dir;
                continue; // No layer folders here, move to next category
            }

            // Now find the REAL target keys inside the layer folder (e.g., eff_eta1_external)
            TIter next_metric_key(layer_blueprint->GetListOfKeys());
            TKey* metric_key = nullptr;

            while ((metric_key = static_cast<TKey*>(next_metric_key()))) {
                std::string metric_name = metric_key->GetName();

                TCanvas* canvas = new TCanvas("c_global", "", 800, 600);
                canvas->cd();

                TMultiGraph* global_multigraph = new TMultiGraph();
                global_multigraph->SetName((metric_name + "_global").c_str());
                global_multigraph->SetTitle((metric_name).c_str());

                bool graph_added = false;

                // STEP 4: Stitch layers together by mapping the path correctly
                for (int layer = 0; layer < 3; ++layer) {
                    std::string scan_folder = "scanned_layer_" + std::to_string(layer);
                    std::string target_path = scan_folder + "/" + analysis_subdir_name + "/layer" + std::to_string(layer) + "/" + metric_name;

                    // Look for either TGraphErrors or a component TGraph
                    TGraph* layer_graph = dynamic_cast<TGraph*>(config_dir->Get(target_path.c_str()));
                    if (!layer_graph) continue;

                    // Clone the graph to avoid ownership issues and add it to the global multigraph
                    TObject* raw_clone = layer_graph->Clone();
                    TGraphErrors* graph_error_clone = dynamic_cast<TGraphErrors*>(raw_clone);
                    if (graph_error_clone) {
                        graph_error_clone->SetName(("graph_layer" + std::to_string(layer)).c_str());
                        global_multigraph->Add(graph_error_clone, "P");
                    } else {
                        TGraph* graph_plain_clone = dynamic_cast<TGraph*>(raw_clone);
                        if (graph_plain_clone) {
                            graph_plain_clone->SetName(("graph_layer" + std::to_string(layer)).c_str());
                            global_multigraph->Add(graph_plain_clone, "P");
                        } else {
                            delete raw_clone;
                        }
                    }
                    graph_added = true;
                }

                // STEP 5: Style and render the output
                if (graph_added) {
                    PlotCategory category = PlotterHelpers::PlotStyler::getPlotCategory(global_multigraph);
                    TClass* cl = TMultiGraph::Class();

                    auto custom_styler = PlotterHelpers::PlotStyler::getCustomStyler(category);

                    if (custom_styler) {
                        custom_styler(global_multigraph, canvas, cl);
                    } else {
                        PlotterHelpers::PlotStyler::styleDefaultPlot(global_multigraph, canvas, cl);
                    }

                    std::filesystem::path export_file = config_output_path / analysis_subdir_name / (metric_name + "_global.pdf");
                    std::filesystem::create_directories(export_file.parent_path());
                    
                    canvas->SaveAs(export_file.string().c_str());
                }

                delete global_multigraph;
                delete canvas;
            }
            delete layer_blueprint;
            delete analysis_dir;
        }
        delete blueprint_dir;
    }

} // namespace BatchExporter
} // namespace PlotterHelpers