#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <initializer_list>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "configParser.hpp"
#include "dataPlotter.hpp"
#include "plotStyler.hpp"
#include "plotBatchExporter.hpp"

#include "utils.hpp"
#include "constants.hpp"

#include "core/hit.hpp"
#include "core/event.hpp"
#include "core/cluster.hpp"
#include "core/track.hpp"
#include "analysis/dataAnalyzer.hpp"

// ==========================================================================================
// Analysis utility/helper namespaces for plotting and calculating statistics
// ==========================================================================================
namespace perFileHelpers {

int remapStrip(int rawStrip) {
    for (const auto& col : columnShifts)
        if (rawStrip >= col.start && rawStrip <= col.end)
            return rawStrip + col.shift;
    return rawStrip;
}

void plotStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* strip_dir = analysis_dir->GetDirectory("strip");
    strip_dir->cd();

    TTree* input_data_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    if (!input_data_tree || input_data_tree->IsZombie()) {
        std::cerr << "Error: Invalid input data tree for analysis." << std::endl;
        return;
    }
    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create strip distributions for each layer and side
    TTreeReader readerInputData(input_data_tree);
    TTreeReader readerProcData(proc_tree);
    TTreeReaderValue<std::vector<int>> raw_time1(readerInputData, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>> raw_time2(readerInputData, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");

    // Create histograms using arrays
    const int nConfigs = 2;
    const char* categories[nConfigs] = {
        "strip_eta1", "strip_eta2"
    };

    std::map<std::string, std::map<int, TH1*>> strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH1F(Form("h1d_%s_layer%d", categories[c], layer),
                            Form("Layer %d;Strip;Hits", layer),
                            24, 0, 24);
            strip_histograms[categories[c]][layer] = hist;
        }
    }

    while (readerInputData.Next() && readerProcData.Next()) {
        for (size_t i = 0; i < strips->size(); ++i) {
            int layer = (*layers)[i];
            int strip = remapStrip((*strips)[i]);
            
            if ((*raw_time1)[i] != 0) strip_histograms["strip_eta1"][layer]->Fill(strip);
            if ((*raw_time2)[i] != 0) strip_histograms["strip_eta2"][layer]->Fill(strip);
        }
    }

    // Write histograms to file and clean up
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            strip_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete strip_histograms[categories[c]][layer];
        }
    }
}

void plotToT(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* tot_dir = analysis_dir->GetDirectory("tot");
    tot_dir->cd();

    TTree* input_data_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    if (!input_data_tree || input_data_tree->IsZombie()) {
        std::cerr << "Error: Invalid input data tree for analysis." << std::endl;
        return;
    }
    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid track reconstruction tree for analysis." << std::endl;
        return;
    }

    TTreeReader readerInputData(input_data_tree);
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> raw_time1(readerInputData, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>> raw_time2(readerInputData, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>> tot1(readerProcData, "proc_tot1");
    TTreeReaderValue<std::vector<int>> tot2(readerProcData, "proc_tot2");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    const int nConfigs = 8;
    const char* categories[nConfigs] = { "tot_eta1", "tot_eta2",
        "tot_eta1_valid1", "tot_eta2_valid1", "tot_eta1_valid2", "tot_eta2_valid2", "tot_eta1_valid_all", "tot_eta2_valid_all" };
    const char* comments[nConfigs] = { "Before Track Reco", "Before Track Reco",
        "After Track Reco (#eta1)", "After Track Reco (#eta2)",
        "After Track Reco (#eta1)", "After Track Reco (#eta2)",
        "After Track Reco (#eta1 & #eta2)", "After Track Reco (#eta1 & #eta2)" };

    const int nBins = 25;
    const float xMin = 0.0;
    const float xMax = 25.0;
    std::map<std::string, std::map<int, TH1*>> strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH1F(Form("h1d_%s_layer%d", categories[c], layer),
                            Form("Layer %d: %s;ToT [ns];Hits", layer, comments[c]),
                            nBins, xMin, xMax);
            strip_histograms[categories[c]][layer] = hist;
        }
    }

    while (readerInputData.Next() && readerProcData.Next() && readerTrackData.Next()) {
        for (size_t i = 0; i < tot1->size(); ++i) {
            int layer = (*layers)[i];

            double tot1_ns = TimeUtils::ticksToTime((*tot1)[i]);
            double tot2_ns = TimeUtils::ticksToTime((*tot2)[i]);
            bool valid_eta1 = (*in_valid_track_eta1)[i];
            bool valid_eta2 = (*in_valid_track_eta2)[i];

            if ((*raw_time1)[i] != 0) strip_histograms["tot_eta1"][layer]->Fill(tot1_ns);
            if ((*raw_time2)[i] != 0) strip_histograms["tot_eta2"][layer]->Fill(tot2_ns);
            if ((*raw_time1)[i] != 0 && valid_eta1) strip_histograms["tot_eta1_valid1"][layer]->Fill(tot1_ns);
            if ((*raw_time2)[i] != 0 && valid_eta1) strip_histograms["tot_eta2_valid1"][layer]->Fill(tot2_ns);
            if ((*raw_time1)[i] != 0 && valid_eta2) strip_histograms["tot_eta1_valid2"][layer]->Fill(tot1_ns);
            if ((*raw_time2)[i] != 0 && valid_eta2) strip_histograms["tot_eta2_valid2"][layer]->Fill(tot2_ns);
            if ((*raw_time1)[i] != 0 && valid_eta1 && valid_eta2) strip_histograms["tot_eta1_valid_all"][layer]->Fill(tot1_ns);
            if ((*raw_time2)[i] != 0 && valid_eta1 && valid_eta2) strip_histograms["tot_eta2_valid_all"][layer]->Fill(tot2_ns);
        }
    }

    struct StackPairing {
        std::string eta1_cat;
        std::string eta2_cat;
        std::string suffix;
        std::string title_modifier;
    };

    std::vector<StackPairing> pairings = {
        {"tot_eta1",           "tot_eta2",           "all",       "Before Track Reco"},
        {"tot_eta1_valid1",    "tot_eta2_valid1",    "valid1",    "After Track Reco (#eta1 Valid)"},
        {"tot_eta1_valid2",    "tot_eta2_valid2",    "valid2",    "After Track Reco (#eta2 Valid)"},
        {"tot_eta1_valid_all", "tot_eta2_valid_all", "valid_all", "After Track Reco (#eta1 and #eta2 Valid)"}
    };

    // Store stacks
    for (int layer : {0, 1, 2}) {
        for (const auto& pair : pairings) {
            auto* h_tot1 = strip_histograms[pair.eta1_cat][layer];
            auto* h_tot2 = strip_histograms[pair.eta2_cat][layer];

            // Style them distinctly inside the stack
            h_tot1->SetLineColor(kBlue);
            h_tot1->SetMarkerColor(kBlue);
            
            h_tot2->SetLineColor(kRed);
            h_tot2->SetMarkerColor(kRed);

            // Construct unique name and descriptive title for each specific stack pairing
            std::string stack_name = Form("h1d_tot_%s_layer%d", pair.suffix.c_str(), layer);
            std::string stack_title = Form("Layer %d: %s", layer, pair.title_modifier.c_str());

            auto* stack = new THStack(stack_name.c_str(), stack_title.c_str());

            // Add to stack (h_tot1 drawn on top of h_tot2 if using "nostack" option later)
            stack->Add(h_tot2);
            stack->Add(h_tot1);

            stack->Write("", TObject::kOverwrite);
            delete stack;
        }
    }

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            strip_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete strip_histograms[categories[c]][layer];
        }
    }
}

void plotDtVsStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* dt_strip_dir = analysis_dir->GetDirectory("dt_strip");
    dt_strip_dir->cd();

    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }
    if (!track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid track reconstruction tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create 2D histograms for dt vs strip for each layer and valid track category
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> dts(readerProcData, "proc_dt_time1_time2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    // Create histograms using arrays
    const int nConfigs = 4;
    const char* suffixes[nConfigs] = {"dt_strip_all", "dt_strip_valid1", "dt_strip_valid2", "dt_strip_valid_all"};
    const char* comments[nConfigs] = {"Before Track Reco", "After Track Reco (#eta1)", "After Track Reco (#eta2)", "After Track Reco (#eta1 and #eta2)"};

    std::map<std::string, std::map<int, TH2*>> dt_strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(Form("h2d_%s_layer%d", suffixes[c], layer),
                            Form("Layer %d: %s;Strip;#Delta#it{t} [Ticks]; Entries", layer, comments[c]),
                            24, 0, 24, 24, -12, 12);
            dt_strip_histograms[suffixes[c]][layer] = hist;
        }
    }

    while (readerProcData.Next() && readerTrackData.Next()) {
        for (size_t i = 0; i < strips->size(); ++i) {
            int layer = (*layers)[i];
            int strip = remapStrip((*strips)[i]);
            double dt_ticks = (*dts)[i];

            dt_strip_histograms["dt_strip_all"][layer]->Fill(strip, dt_ticks);
            if ((*in_valid_track_eta1)[i]) {
                dt_strip_histograms["dt_strip_valid1"][layer]->Fill(strip, dt_ticks);
            }
            if ((*in_valid_track_eta2)[i]) {
                dt_strip_histograms["dt_strip_valid2"][layer]->Fill(strip, dt_ticks);
            }
            if ((*in_valid_track_eta1)[i] && (*in_valid_track_eta2)[i]) {
                dt_strip_histograms["dt_strip_valid_all"][layer]->Fill(strip, dt_ticks);
            }
        }
    }

    // Write histograms to file and clean up
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            dt_strip_histograms[suffixes[c]][layer]->Write("", TObject::kOverwrite);
            delete dt_strip_histograms[suffixes[c]][layer];
        }
    }
}

void plotToTVsStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* tot_strip_dir = analysis_dir->GetDirectory("tot_strip");
    tot_strip_dir->cd();

    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }
    if (!track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid track reconstruction tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create 2D histograms for dt vs strip for each layer and valid track category
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> proc_tot1(readerProcData, "proc_tot1");
    TTreeReaderValue<std::vector<int>> proc_tot2(readerProcData, "proc_tot2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    // Create histograms using arrays
    const int nConfigs = 8;
    const char* suffixes[nConfigs] = {"tot1_strip_all", "tot2_strip_all", "tot1_strip_valid1", "tot2_strip_valid1",
        "tot1_strip_valid2", "tot2_strip_valid2", "tot1_strip_valid_all", "tot2_strip_valid_all"};
    const char* comments[nConfigs] = {"Before Track Reco", "Before Track Reco", "After Track Reco (#eta1)", "After Track Reco (#eta1)",
        "After Track Reco (#eta2)", "After Track Reco (#eta2)", "After Track Reco (#eta1 and #eta2)", "After Track Reco (#eta1 and #eta2)"};

    std::map<std::string, std::map<int, TH2*>> tot_strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(Form("h2d_%s_layer%d", suffixes[c], layer),
                            Form("Layer %d: %s;Strip;ToT [ns]; Entries", layer, comments[c]),
                            24, 0, 24, 34, 1, 35);
            tot_strip_histograms[suffixes[c]][layer] = hist;
        }
    }

    while (readerProcData.Next() && readerTrackData.Next()) {
        for (size_t i = 0; i < strips->size(); ++i) {
            int layer = (*layers)[i];
            
            int strip = remapStrip((*strips)[i]);
            int tot1 = (*proc_tot1)[i];
            int tot2 = (*proc_tot2)[i];

            
            tot_strip_histograms["tot1_strip_all"][layer]->Fill(strip, tot1);
            tot_strip_histograms["tot2_strip_all"][layer]->Fill(strip, tot2);
            if ((*in_valid_track_eta1)[i]) {
                tot_strip_histograms["tot1_strip_valid1"][layer]->Fill(strip, tot1);
                tot_strip_histograms["tot2_strip_valid1"][layer]->Fill(strip, tot2);
            }
            if ((*in_valid_track_eta2)[i]) {
                tot_strip_histograms["tot1_strip_valid2"][layer]->Fill(strip, tot1);
                tot_strip_histograms["tot2_strip_valid2"][layer]->Fill(strip, tot2);
            }
            if ((*in_valid_track_eta1)[i] && (*in_valid_track_eta2)[i]) {
                tot_strip_histograms["tot1_strip_valid_all"][layer]->Fill(strip, tot1);
                tot_strip_histograms["tot2_strip_valid_all"][layer]->Fill(strip, tot2);
            }
        }
    }

    // Write histograms to file and clean up
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {   
            tot_strip_histograms[suffixes[c]][layer]->Write("", TObject::kOverwrite);
            delete tot_strip_histograms[suffixes[c]][layer];
        }
    }

}

void plotMultiplicityAndDelayVsStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* mult_strip_dir = analysis_dir->GetDirectory("multiplicity_strip");
    TDirectory* delay_strip_dir = analysis_dir->GetDirectory("delay_strip");

    TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!proc_tree || proc_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed data tree for analysis." << std::endl;
        return;
    }
    if (!track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid track reconstruction tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create 2D histograms for dt vs strip for each layer and valid track category
    TTreeReader readerInputData(input_tree);
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> raw_time1(readerInputData, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>> raw_time2(readerInputData, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>> rise(readerInputData, "hit_rise");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> time1(readerProcData, "proc_time1");
    TTreeReaderValue<std::vector<int>> time2(readerProcData, "proc_time2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    // Create histograms using arrays
    const int nConfigs = 8;
    const char* categories[nConfigs] = {
        "mult1_all", "mult2_all",           // Multiplicity for all hits based on time1 and time2
        "mult1_valid1", "mult2_valid1",     // Multiplicity for hits in valid tracks on eta1 side based on time1 and time2
        "mult1_valid2", "mult2_valid2",      // Multiplicity for hits in valid tracks on eta2 side based on time1 and time2
        "mult1_valid_all", "mult2_valid_all" // Multiplicity for hits in valid tracks on both sides based on time1 and time2
    };
    const char* comments[nConfigs] = {
        "Before Track Reco", "Before Track Reco",
        "After Track Reco (#eta1)", "After Track Reco (#eta1)",
        "After Track Reco (#eta2)", "After Track Reco (#eta2)",
        "After Track Reco (#eta1 and #eta2)", "After Track Reco (#eta1 and #eta2)"
    };

    std::map<std::string, std::map<int, TH2*>> multiplicity_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(Form("h2d_%s_layer%d", categories[c], layer),
                            Form("Layer %d: %s;Strip;Multiplicity; Entries", layer, comments[c]),
                            24, 0, 24, 9, 1, 10);
            multiplicity_histograms[categories[c]][layer] = hist;
        }
    }
    std::map<std::string, std::map<int, TH2*>> delay_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(Form("h2d_delay_%s_layer%d", categories[c], layer),
                            Form("Layer %d: %s;Strip;Delay from First Hit [Ticks]; Entries", layer, comments[c]),
                            24, 0, 24, 100, 0, 100);
            delay_histograms[categories[c]][layer] = hist;
        }
    }

    while (readerInputData.Next() && readerProcData.Next() && readerTrackData.Next()) {

        // Create lookup tables for each entry
        std::map<int, std::map<int, int>> counts[8];    // layer -> strip -> count (for each of the 8 categories)
        std::map<int, std::map<int, std::vector<int>>> delays[8];

        for (size_t i = 0; i < strips->size(); ++i) {

            if ((*rise)[i] == 0) continue; // Skip falling edge hits

            int layer = (*layers)[i];
            int strip = remapStrip((*strips)[i]);

            if ((*raw_time1)[i] != 0) {
                counts[0][layer][strip]++;
                delays[0][layer][strip].push_back((*time1)[i]);
                if ((*in_valid_track_eta1)[i]) {
                    counts[2][layer][strip]++;
                    delays[2][layer][strip].push_back((*time1)[i]);
                }
                if ((*in_valid_track_eta2)[i]) {
                    counts[4][layer][strip]++;
                    delays[4][layer][strip].push_back((*time1)[i]);
                }
                if ((*in_valid_track_eta1)[i] && (*in_valid_track_eta2)[i]) {
                    counts[6][layer][strip]++;
                    delays[6][layer][strip].push_back((*time1)[i]);
                }
            }
            if ((*raw_time2)[i] != 0) {
                counts[1][layer][strip]++;
                delays[1][layer][strip].push_back((*time2)[i]);
                if ((*in_valid_track_eta1)[i]) {
                    counts[3][layer][strip]++;
                    delays[3][layer][strip].push_back((*time2)[i]);
                }
                if ((*in_valid_track_eta2)[i]) {
                    counts[5][layer][strip]++;
                    delays[5][layer][strip].push_back((*time2)[i]);
                }
                if ((*in_valid_track_eta1)[i] && (*in_valid_track_eta2)[i]) {
                    counts[7][layer][strip]++;
                    delays[7][layer][strip].push_back((*time2)[i]);
                }
            }
        }

        // Fill histograms for this entry based on the counts and delays in the lookup tables
        for (int c = 0; c < nConfigs; ++c) {
            for (int layer : {0, 1, 2}) {
                for (const auto& [strip, count] : counts[c][layer]) {
                    multiplicity_histograms[categories[c]][layer]->Fill(strip, count);
                }

                for (const auto& [strip, delay_vec] : delays[c][layer]) {
                    if (!delay_vec.empty()) {
                        auto min_delay = std::min_element(delay_vec.begin(), delay_vec.end());
                        size_t min_index = std::distance(delay_vec.begin(), min_delay);
                        for (size_t i = 0; i < delay_vec.size(); ++i) {
                            if (i == min_index) continue; // Skip the first hit (minimum delay)
                            delay_histograms[categories[c]][layer]->Fill(strip, delay_vec[i] - *min_delay);
                        }
                    }
                }
            }
        }
    }

    // Write histograms to file and clean up
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            mult_strip_dir->cd();
            multiplicity_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete multiplicity_histograms[categories[c]][layer];
            delay_strip_dir->cd();
            delay_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete delay_histograms[categories[c]][layer];
        }
    }
}

const std::array<std::string, LAYER_PAIR_COUNT> pair_labels = {"Layer 0-1", "Layer 0-2", "Layer 1-2"};

void plotToFs(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    if (!analysis_dir) return;

    TDirectory* tof_dir = analysis_dir->GetDirectory("tof");
    if (!tof_dir) return;
    tof_dir->cd();

    TTree* track_tree = input_file->Get<TTree>("TrackReconstruction");
    if (!track_tree) {
        std::cerr << "Error: TrackReconstruction tree not found in file." << std::endl;
        return;
    }

    TTreeReader readerTrackData(track_tree);

    std::array<std::unique_ptr<TTreeReaderValue<std::vector<int>>>, LAYER_PAIR_COUNT> tof1_readers;
    std::array<std::unique_ptr<TTreeReaderValue<std::vector<int>>>, LAYER_PAIR_COUNT> tof2_readers;

    std::array<TH1F*, LAYER_PAIR_COUNT> h_tof1;
    std::array<TH1F*, LAYER_PAIR_COUNT> h_tof2;

    for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
        std::string branch1 = "track_time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta1";
        std::string branch2 = "track_time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta2";

        tof1_readers[i] = std::make_unique<TTreeReaderValue<std::vector<int>>>(readerTrackData, branch1.c_str());
        tof2_readers[i] = std::make_unique<TTreeReaderValue<std::vector<int>>>(readerTrackData, branch2.c_str());

        h_tof1[i] = new TH1F(Form("h1d_tof_layer_%s_eta1", LAYER_PAIR_SUFFIXES[i].c_str()),
                             Form("Side #eta1: %s;ToF [Ticks];Entries", pair_labels[i].c_str()), 16, -8, 8);
        h_tof2[i] = new TH1F(Form("h1d_tof_layer_%s_eta2", LAYER_PAIR_SUFFIXES[i].c_str()),
                             Form("Side #eta2: %s;ToF [Ticks];Entries", pair_labels[i].c_str()), 16, -8, 8);
    }

    while (readerTrackData.Next()) {
        for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {

            if (tof1_readers[i]->GetSetupStatus() == 0) {
                for (int t : **tof1_readers[i]) {
                    h_tof1[i]->Fill(t);
                }
            }

            if (tof2_readers[i]->GetSetupStatus() == 0) {
                for (int t : **tof2_readers[i]) {
                    h_tof2[i]->Fill(t);
                }
            }
        }
    }

    for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
        h_tof1[i]->Write("", TObject::kOverwrite);
        h_tof2[i]->Write("", TObject::kOverwrite);

        delete h_tof1[i];
        delete h_tof2[i];
    }
}

}

namespace summaryHelpers {

void setupBranches(TTree* summary_tree, MeasurementMetadata& metadata, MeasurementData& data) {
    summary_tree->Branch("name", &metadata.name);
    summary_tree->Branch("measurement_type", &metadata.measurement_type);
    summary_tree->Branch("group_name", &metadata.group_name);
    summary_tree->Branch("mixture", &metadata.mixture);
    summary_tree->Branch("source", &metadata.source);
    summary_tree->Branch("filter", &metadata.filter);
    summary_tree->Branch("lv_setting", &metadata.lv_setting);
    summary_tree->Branch("scanned_layer", &metadata.scanned_layer);
    summary_tree->Branch("scanned_hv", &metadata.scanned_hv);
    summary_tree->Branch("other_hv", &metadata.other_hv);

    summary_tree->Branch("eff_eta1_external", &data.efficiency_results.eta1_efficiency_external, "eff_eta1_external[3]/D");
    summary_tree->Branch("eff_eta2_external", &data.efficiency_results.eta2_efficiency_external, "eff_eta2_external[3]/D");
    summary_tree->Branch("eff_or_external", &data.efficiency_results.eta_or_efficiency_external, "eff_or_external[3]/D");
    summary_tree->Branch("eff_and_external", &data.efficiency_results.eta_and_efficiency_external, "eff_and_external[3]/D");
    summary_tree->Branch("eff_eta1_external_error", &data.efficiency_results.eta1_efficiency_external_error, "eff_eta1_external_error[6]/D");
    summary_tree->Branch("eff_eta2_external_error", &data.efficiency_results.eta2_efficiency_external_error, "eff_eta2_external_error[6]/D");
    summary_tree->Branch("eff_or_external_error", &data.efficiency_results.eta_or_efficiency_external_error, "eff_or_external_error[6]/D");
    summary_tree->Branch("eff_and_external_error", &data.efficiency_results.eta_and_efficiency_external_error, "eff_and_external_error[6]/D");

    summary_tree->Branch("eff_eta1_rpc", &data.efficiency_results.eta1_efficiency_rpc, "eff_eta1_rpc[3]/D");
    summary_tree->Branch("eff_eta2_rpc", &data.efficiency_results.eta2_efficiency_rpc, "eff_eta2_rpc[3]/D");
    summary_tree->Branch("eff_or_rpc", &data.efficiency_results.eta_or_efficiency_rpc, "eff_or_rpc[3]/D");
    summary_tree->Branch("eff_and_rpc", &data.efficiency_results.eta_and_efficiency_rpc, "eff_and_rpc[3]/D");
    summary_tree->Branch("eff_eta1_rpc_error", &data.efficiency_results.eta1_efficiency_rpc_error, "eff_eta1_rpc_error[6]/D");
    summary_tree->Branch("eff_eta2_rpc_error", &data.efficiency_results.eta2_efficiency_rpc_error, "eff_eta2_rpc_error[6]/D");
    summary_tree->Branch("eff_or_rpc_error", &data.efficiency_results.eta_or_efficiency_rpc_error, "eff_or_rpc_error[6]/D");
    summary_tree->Branch("eff_and_rpc_error", &data.efficiency_results.eta_and_efficiency_rpc_error, "eff_and_rpc_error[6]/D");

    summary_tree->Branch("track_eff_eta1_external", &data.efficiency_results_tracks.eta1_efficiency_external, "track_eff_eta1_external[3]/D");
    summary_tree->Branch("track_eff_eta2_external", &data.efficiency_results_tracks.eta2_efficiency_external, "track_eff_eta2_external[3]/D");
    summary_tree->Branch("track_eff_or_external", &data.efficiency_results_tracks.eta_or_efficiency_external, "track_eff_or_external[3]/D");
    summary_tree->Branch("track_eff_and_external", &data.efficiency_results_tracks.eta_and_efficiency_external, "track_eff_and_external[3]/D");
    summary_tree->Branch("track_eff_eta1_external_error", &data.efficiency_results_tracks.eta1_efficiency_external_error, "track_eff_eta1_external_error[6]/D");
    summary_tree->Branch("track_eff_eta2_external_error", &data.efficiency_results_tracks.eta2_efficiency_external_error, "track_eff_eta2_external_error[6]/D");
    summary_tree->Branch("track_eff_or_external_error", &data.efficiency_results_tracks.eta_or_efficiency_external_error, "track_eff_or_external_error[6]/D");
    summary_tree->Branch("track_eff_and_external_error", &data.efficiency_results_tracks.eta_and_efficiency_external_error, "track_eff_and_external_error[6]/D");

    summary_tree->Branch("track_eff_eta1_rpc", &data.efficiency_results_tracks.eta1_efficiency_rpc, "track_eff_eta1_rpc[3]/D");
    summary_tree->Branch("track_eff_eta2_rpc", &data.efficiency_results_tracks.eta2_efficiency_rpc, "track_eff_eta2_rpc[3]/D");
    summary_tree->Branch("track_eff_or_rpc", &data.efficiency_results_tracks.eta_or_efficiency_rpc, "track_eff_or_rpc[3]/D");
    summary_tree->Branch("track_eff_and_rpc", &data.efficiency_results_tracks.eta_and_efficiency_rpc, "track_eff_and_rpc[3]/D");
    summary_tree->Branch("track_eff_eta1_rpc_error", &data.efficiency_results_tracks.eta1_efficiency_rpc_error, "track_eff_eta1_rpc_error[6]/D");
    summary_tree->Branch("track_eff_eta2_rpc_error", &data.efficiency_results_tracks.eta2_efficiency_rpc_error, "track_eff_eta2_rpc_error[6]/D");
    summary_tree->Branch("track_eff_or_rpc_error", &data.efficiency_results_tracks.eta_or_efficiency_rpc_error, "track_eff_or_rpc_error[6]/D");
    summary_tree->Branch("track_eff_and_rpc_error", &data.efficiency_results_tracks.eta_and_efficiency_rpc_error, "track_eff_and_rpc_error[6]/D");

    summary_tree->Branch("avg_cluster_size_eta1", &data.cluster_size_results.avg_cluster_size_eta1);
    summary_tree->Branch("avg_cluster_size_eta2", &data.cluster_size_results.avg_cluster_size_eta2);
    summary_tree->Branch("avg_cluster_size_eta1_layers", &data.cluster_size_results.avg_cluster_size_eta1_layers, "avg_cluster_size_eta1_layers[3]/D");
    summary_tree->Branch("avg_cluster_size_eta2_layers", &data.cluster_size_results.avg_cluster_size_eta2_layers, "avg_cluster_size_eta2_layers[3]/D");
    summary_tree->Branch("avg_cluster_size_eta1_error", &data.cluster_size_results.avg_cluster_size_eta1_error, "avg_cluster_size_eta1_error[2]/D");
    summary_tree->Branch("avg_cluster_size_eta2_error", &data.cluster_size_results.avg_cluster_size_eta2_error, "avg_cluster_size_eta2_error[2]/D");
    summary_tree->Branch("avg_cluster_size_eta1_layers_error", &data.cluster_size_results.avg_cluster_size_eta1_layers_error, "avg_cluster_size_eta1_layers_error[6]/D");
    summary_tree->Branch("avg_cluster_size_eta2_layers_error", &data.cluster_size_results.avg_cluster_size_eta2_layers_error, "avg_cluster_size_eta2_layers_error[6]/D");

    summary_tree->Branch("noise_rate", &data.noise_rate_results.noise_rate);
    summary_tree->Branch("noise_rate_error", &data.noise_rate_results.noise_rate_error, Form("noise_rate_error[%d]/D", 2));

    summary_tree->Branch("noise_rate_eta1", &data.noise_rate_results.noise_rate_eta1, Form("noise_rate_eta1[%d]/D", 3));
    summary_tree->Branch("noise_rate_eta2", &data.noise_rate_results.noise_rate_eta2, Form("noise_rate_eta2[%d]/D", 3));
    summary_tree->Branch("noise_rate_eta1_error", &data.noise_rate_results.noise_rate_eta1_error, Form("noise_rate_eta1_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("noise_rate_eta2_error", &data.noise_rate_results.noise_rate_eta2_error, Form("noise_rate_eta2_error[%d]/D", LAYER_COUNT * 2));

    summary_tree->Branch("noise_rate_strips_eta1", &data.noise_rate_results.noise_rate_strips_eta1, Form("noise_rate_strips_eta1[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("noise_rate_strips_eta2", &data.noise_rate_results.noise_rate_strips_eta2, Form("noise_rate_strips_eta2[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("noise_rate_strips_eta1_error", &data.noise_rate_results.noise_rate_strips_eta1_error, Form("noise_rate_strips_eta1_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("noise_rate_strips_eta2_error", &data.noise_rate_results.noise_rate_strips_eta2_error, Form("noise_rate_strips_eta2_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("avg_tot_eta1", &data.tot_results.avg_tot_eta1, Form("avg_tot_eta1[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_tot_eta2", &data.tot_results.avg_tot_eta2, Form("avg_tot_eta2[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_tot_eta1_error", &data.tot_results.avg_tot_eta1_error, Form("avg_tot_eta1_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("avg_tot_eta2_error", &data.tot_results.avg_tot_eta2_error, Form("avg_tot_eta2_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("track_avg_tot_eta1", &data.tot_results_tracks.avg_tot_eta1, Form("track_avg_tot_eta1[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("track_avg_tot_eta2", &data.tot_results_tracks.avg_tot_eta2, Form("track_avg_tot_eta2[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("track_avg_tot_eta1_error", &data.tot_results_tracks.avg_tot_eta1_error, Form("track_avg_tot_eta1_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("track_avg_tot_eta2_error", &data.tot_results_tracks.avg_tot_eta2_error, Form("track_avg_tot_eta2_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("avg_multiplicity_eta1", &data.multiplicity_results.avg_multiplicity_eta1, Form("avg_multiplicity_eta1[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_multiplicity_eta2", &data.multiplicity_results.avg_multiplicity_eta2, Form("avg_multiplicity_eta2[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_multiplicity_eta1_error", &data.multiplicity_results.avg_multiplicity_eta1_error, Form("avg_multiplicity_eta1_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("avg_multiplicity_eta2_error", &data.multiplicity_results.avg_multiplicity_eta2_error, Form("avg_multiplicity_eta2_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("track_avg_multiplicity_eta1", &data.multiplicity_results_tracks.avg_multiplicity_eta1, Form("track_avg_multiplicity_eta1[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("track_avg_multiplicity_eta2", &data.multiplicity_results_tracks.avg_multiplicity_eta2, Form("track_avg_multiplicity_eta2[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("track_avg_multiplicity_eta1_error", &data.multiplicity_results_tracks.avg_multiplicity_eta1_error, Form("track_avg_multiplicity_eta1_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("track_avg_multiplicity_eta2_error", &data.multiplicity_results_tracks.avg_multiplicity_eta2_error, Form("track_avg_multiplicity_eta2_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("time_of_flight_eta1", &data.tof_results.time_of_flight_eta1);
    summary_tree->Branch("time_of_flight_eta2", &data.tof_results.time_of_flight_eta2);

    summary_tree->Branch("avg_time_of_flight_eta1", &data.tof_results.avg_time_of_flight_eta1);
    summary_tree->Branch("avg_time_of_flight_eta2", &data.tof_results.avg_time_of_flight_eta2);
    summary_tree->Branch("avg_time_of_flight_eta1_error", &data.tof_results.avg_time_of_flight_eta1_error, Form("avg_time_of_flight_eta1_error[%d]/D", 2));
    summary_tree->Branch("avg_time_of_flight_eta2_error", &data.tof_results.avg_time_of_flight_eta2_error, Form("avg_time_of_flight_eta2_error[%d]/D", 2));

    summary_tree->Branch("time_resolution_eta1", &data.time_resolution_results.time_resolution_eta1);
    summary_tree->Branch("time_resolution_eta2", &data.time_resolution_results.time_resolution_eta2);
    summary_tree->Branch("time_resolution_eta1_error", &data.time_resolution_results.time_resolution_eta1_error, Form("time_resolution_eta1_error[%d]/D", 2));
    summary_tree->Branch("time_resolution_eta2_error", &data.time_resolution_results.time_resolution_eta2_error, Form("time_resolution_eta2_error[%d]/D", 2));
}

inline void requireEqualSizes(std::initializer_list<std::pair<std::string, size_t>> named_sizes) {
    if (named_sizes.size() <= 1) return;

    auto it = named_sizes.begin();
    const std::string& expected_name = it->first;
    size_t expected_size = it->second;

    for (++it; it != named_sizes.end(); ++it) {
        if (it->second != expected_size) {
            throw std::runtime_error(
                "Data size mismatch detected! '" + expected_name +
                "' has size " + std::to_string(expected_size) +
                ", but '" + it->first +
                "' has size " + std::to_string(it->second) + "."
            );
        }
    }
}

void getRate(TFile* input_file, NoiseRateResults& rate_results) {
    
    auto calculateRate = [](long long hits, long long events, int n_strips, int n_layers)
        -> std::pair<double, double> {

        if (events <= 0 || n_strips <= 0 || n_layers <= 0) {
            return {0.0, 0.0};
        }

        // Total time-area exposure = events * window * sensitive area
        const double total_exposure = events * TRIGGER_TIME_WINDOW * n_layers * n_strips * STRIP_WIDTH_CM * DETECTOR_LENGTH_CM;

        const double rate = static_cast<double>(hits) / total_exposure;
        const double rate_error = (hits > 0) ? (std::sqrt(static_cast<double>(hits)) / total_exposure) : 0.0;

        return {rate, rate_error};
    };

    TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    TTree* processed_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    if (!input_tree || !processed_tree) {
        std::cerr << "Error: One or more required trees not found in file " << input_file->GetName() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    TTreeReader readerInputData(input_tree);
    TTreeReader readerProcData(processed_tree);

    TTreeReaderValue<int> hits(readerProcData, "n_hits");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> time1(readerProcData, "proc_time1");
    TTreeReaderValue<std::vector<int>> time2(readerProcData, "proc_time2");

    // Global counters
    long long total_hits = 0;
    long long event_count = 0;
    std::unordered_set<int> unique_strips;
    std::unordered_set<int> unique_layers;

    // Eta/Layer-specific counters
    std::array<long long, 3> total_hits_eta1 = {0, 0, 0};
    std::array<long long, 3> total_hits_eta2 = {0, 0, 0};
    std::array<std::unordered_set<int>, 3> unique_strips_eta1;
    std::array<std::unordered_set<int>, 3> unique_strips_eta2;

    // Eta/Layer/Strip-specific counters
    long long hits_per_strip_eta1[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};
    long long hits_per_strip_eta2[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};

    while (readerInputData.Next() && readerProcData.Next()) {
        event_count += 1;

        // Process Global Rates
        if (hits.GetSetupStatus() == 0) {
            total_hits += static_cast<long long>(*hits);
        }
        if (strips.GetSetupStatus() == 0 && layers.GetSetupStatus() == 0) {
            for (int strip : *strips) {
                unique_strips.insert(strip);
            }
            for (int layer : *layers) {
                unique_layers.insert(layer);
            }
        }

        // Process Eta/Layer-specific Rates
        if (time1.GetSetupStatus() == 0 && time2.GetSetupStatus() == 0) {
            const size_t n_hits = std::min({time1->size(), time2->size()});
            for (size_t i = 0; i < n_hits; ++i) {
                int layer = (*layers)[i];
                int strip = perFileHelpers::remapStrip((*strips)[i]);

                if ((*time1)[i] != 0) {
                    total_hits_eta1[layer] += 1;
                    hits_per_strip_eta1[layer][strip]++;
                    unique_strips_eta1[layer].insert(strip);
                }
                if ((*time2)[i] != 0) {
                    total_hits_eta2[layer] += 1;
                    hits_per_strip_eta2[layer][strip]++;
                    unique_strips_eta2[layer].insert(strip);
                }
            }
        }
    }

    // Calculate final global rates
    int n_strips = static_cast<int>(unique_strips.size());
    int n_layers = static_cast<int>(unique_layers.size());

    std::tie(rate_results.noise_rate, rate_results.noise_rate_error) = 
            calculateRate(total_hits, event_count, n_strips, n_layers);

    // Calculate final Eta/Layer-specific rates
    for (int layer = 0; layer < LAYER_COUNT; ++layer) {
        const int n_strips_eta1 = static_cast<int>(unique_strips_eta1[layer].size());
        const int n_strips_eta2 = static_cast<int>(unique_strips_eta2[layer].size());

        std::tie(rate_results.noise_rate_eta1[layer], rate_results.noise_rate_eta1_error[layer]) = 
            calculateRate(total_hits_eta1[layer], event_count, n_strips_eta1, 1);

        std::tie(rate_results.noise_rate_eta2[layer], rate_results.noise_rate_eta2_error[layer]) = 
            calculateRate(total_hits_eta2[layer], event_count, n_strips_eta2, 1);

        // Calculate final Strip-specific rates
        for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {
            std::tie(rate_results.noise_rate_strips_eta1[layer][strip], 
                     rate_results.noise_rate_strips_eta1_error[layer][strip]) = 
                calculateRate(hits_per_strip_eta1[layer][strip], event_count, 1, 1);

            std::tie(rate_results.noise_rate_strips_eta2[layer][strip], 
                     rate_results.noise_rate_strips_eta2_error[layer][strip]) = 
                calculateRate(hits_per_strip_eta2[layer][strip], event_count, 1, 1);
        }
    }
}

void getAverageToT(TFile* input_file, ToTResults& tot_results, bool in_valid_track_only) {
    TTree* processed_tree = input_file->Get<TTree>("ProcessedData");
    TTree* track_tree = input_file->Get<TTree>("TrackReconstruction");

    if (!processed_tree || !track_tree) {
        std::cerr << "Error: Required trees not found in file " << input_file->GetName() << "\n";
        std::exit(EXIT_FAILURE);
    }

    TTreeReader reader_processed(processed_tree);
    TTreeReader reader_track(track_tree);

    TTreeReaderValue<std::vector<int>>  tot1(reader_processed, "proc_tot1");
    TTreeReaderValue<std::vector<int>>  tot2(reader_processed, "proc_tot2");
    TTreeReaderValue<std::vector<int>>  strips(reader_processed, "proc_strip");
    TTreeReaderValue<std::vector<int>>  layers(reader_processed, "proc_layer");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(reader_track, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(reader_track, "in_valid_track_eta2");

    // Local accumulators
    Accumulator eta1[LAYER_COUNT][STRIPS_PER_LAYER];
    Accumulator eta2[LAYER_COUNT][STRIPS_PER_LAYER];

    while (reader_processed.Next() && reader_track.Next()) {
        if (tot1.GetSetupStatus() != 0 || in_valid_track_eta1.GetSetupStatus() != 0) continue;

        /* WIP: Figuring out why track reco leads to fewer entries per event than the processed data, which is causing size mismatches
        requireEqualSizes({
            {"tot1", tot1->size()},
            {"tot2", tot2->size()},
            {"strips", strips->size()},
            {"layers", layers->size()},
            {"in_valid_track_eta1", in_valid_track_eta1->size()},
            {"in_valid_track_eta2", in_valid_track_eta2->size()}
        });
        */

        const size_t n_hits = std::min({tot1->size(), tot2->size(), strips->size(), layers->size(),
            in_valid_track_eta1->size(), in_valid_track_eta2->size()});

        for (size_t i = 0; i < n_hits; ++i) {
            int layer = (*layers)[i];
            int strip = perFileHelpers::remapStrip((*strips)[i]);

            if (layer < 0 || layer >= LAYER_COUNT || strip < 0 || strip >= STRIPS_PER_LAYER) continue;

            // Accumulate for all hits or only for hits that are part of valid tracks based on the flag
            if (!in_valid_track_only) {
                if ((*tot1)[i] > 0) {
                    eta1[layer][strip].sum += (*tot1)[i];
                    eta1[layer][strip].hits++;
                }
                if ((*tot2)[i] > 0) {
                    eta2[layer][strip].sum += (*tot2)[i];
                    eta2[layer][strip].hits++;
                }
            } else {
                if ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]) {
                    if ((*tot1)[i] > 0) {
                        eta1[layer][strip].sum += (*tot1)[i];
                        eta1[layer][strip].hits++;
                    }
                    if ((*tot2)[i] > 0) {
                        eta2[layer][strip].sum += (*tot2)[i];
                        eta2[layer][strip].hits++;
                    }
                }
            }
        }
    }

    auto assignAverageAndError = [](const Accumulator& acc, double& avg, ErrorRange& err) {
        if (acc.hits == 0) {
            avg = 0.0;
            err = ErrorRange(0.0);
        } else {
            avg = static_cast<double>(acc.sum) / acc.hits;
            err = ErrorRange(std::sqrt(avg) / acc.hits);
        }
    };

    for (int layer = 0; layer < LAYER_COUNT; ++layer) {
        for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {

            assignAverageAndError(eta1[layer][strip],
                                  tot_results.avg_tot_eta1[layer][strip],
                                  tot_results.avg_tot_eta1_error[layer][strip]);

            assignAverageAndError(eta2[layer][strip],
                                  tot_results.avg_tot_eta2[layer][strip],
                                  tot_results.avg_tot_eta2_error[layer][strip]);
        }
    }
}

void getAverageMultiplicity(TFile* input_file, MultiplicityResults& mult_results, bool in_valid_track_only) {
    TTree* input_tree = input_file->Get<TTree>("InputData");
    TTree* processed_tree = input_file->Get<TTree>("ProcessedData");
    TTree* track_tree = input_file->Get<TTree>("TrackReconstruction");

    if (!input_tree || !processed_tree || !track_tree) {
        std::cerr << "Error: Required trees not found in file " << input_file->GetName() << "\n";
        std::exit(EXIT_FAILURE);
    }

    TTreeReader reader_input(input_tree);
    TTreeReader reader_processed(processed_tree);
    TTreeReader reader_track(track_tree);

    TTreeReaderValue<std::vector<int>>  raw_time1(reader_input, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>>  raw_time2(reader_input, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>>  rise(reader_input, "hit_rise");
    TTreeReaderValue<std::vector<int>>  strips(reader_processed, "proc_strip");
    TTreeReaderValue<std::vector<int>>  layers(reader_processed, "proc_layer");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(reader_track, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(reader_track, "in_valid_track_eta2");

    // Local counters
    int n_events = 0;
    int eta1_per_layer_strip[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};
    int eta2_per_layer_strip[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};

    // Track how many events actually had at least 1 hit on a specific strip
    int active_events_eta1[LAYER_COUNT][STRIPS_PER_LAYER] = {0};
    int active_events_eta2[LAYER_COUNT][STRIPS_PER_LAYER] = {0};

    while (reader_input.Next() && reader_processed.Next() && reader_track.Next()) {
        if (in_valid_track_eta1.GetSetupStatus() != 0 || in_valid_track_eta2.GetSetupStatus() != 0) continue;

        const size_t n_hits = std::min({strips->size(), layers->size(), in_valid_track_eta1->size(), in_valid_track_eta2->size()});

        // Temporary flags to see if a strip fired in this event
        bool hit_in_event_eta1[LAYER_COUNT][STRIPS_PER_LAYER] = {false};
        bool hit_in_event_eta2[LAYER_COUNT][STRIPS_PER_LAYER] = {false};

        for (size_t i = 0; i < n_hits; ++i) {
            if ((*rise)[i] == 0) continue; // Skip falling edge hits

            int layer = (*layers)[i];
            int strip = perFileHelpers::remapStrip((*strips)[i]);

            if (layer < 0 || layer >= LAYER_COUNT || strip < 0 || strip >= STRIPS_PER_LAYER) continue;

            // Accumulate for all hits or only for hits that are part of valid tracks based on the flag
            if (!in_valid_track_only) {
                if ((*raw_time1)[i] != 0) {
                    eta1_per_layer_strip[layer][strip]++;
                    hit_in_event_eta1[layer][strip] = true;
                }
                if ((*raw_time2)[i] != 0) {
                    eta2_per_layer_strip[layer][strip]++;
                    hit_in_event_eta2[layer][strip] = true;
                }
            } else {
                if ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]) {
                    if ((*raw_time1)[i] != 0) {
                        eta1_per_layer_strip[layer][strip]++;
                        hit_in_event_eta1[layer][strip] = true;
                    }
                    if ((*raw_time2)[i] != 0) {
                        eta2_per_layer_strip[layer][strip]++;
                        hit_in_event_eta2[layer][strip] = true;
                    }
                }
            }
        }
        for (int layer = 0; layer < LAYER_COUNT; ++layer) {
            for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {
                if (hit_in_event_eta1[layer][strip]) active_events_eta1[layer][strip]++;
                if (hit_in_event_eta2[layer][strip]) active_events_eta2[layer][strip]++;
            }
        }
        n_events++;
    }

    auto assignAverageAndError = [](int count, int active_events, double& avg, ErrorRange& err) {
        if (active_events == 0) {
            avg = 0.0;
            err = ErrorRange(0.0);
        } else {
            avg = static_cast<double>(count) / active_events;
            err = ErrorRange(std::sqrt(count) / active_events);
        }
    };

    for (int layer = 0; layer < LAYER_COUNT; ++layer) {
        for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {

            assignAverageAndError(eta1_per_layer_strip[layer][strip],
                                  active_events_eta1[layer][strip],
                                  mult_results.avg_multiplicity_eta1[layer][strip],
                                  mult_results.avg_multiplicity_eta1_error[layer][strip]);

            assignAverageAndError(eta2_per_layer_strip[layer][strip],
                                  active_events_eta2[layer][strip],
                                  mult_results.avg_multiplicity_eta2[layer][strip],
                                  mult_results.avg_multiplicity_eta2_error[layer][strip]);
        }
    }
}

void processToF(TFile* input_file, ToFResults& tof_results, TimeResolutionResults& time_resolution_results) {
    TTree* track_tree = input_file->Get<TTree>("TrackReconstruction");

    if (!track_tree) {
        std::cerr << "Error: Required trees not found in file " << input_file->GetName() << "\n";
        std::exit(EXIT_FAILURE);
    }

    TTreeReader reader_track(track_tree);
    TTreeReaderValue<std::vector<int>> tof1(reader_track, "track_time_of_flight_eta1");
    TTreeReaderValue<std::vector<int>> tof2(reader_track, "track_time_of_flight_eta2");

    TH1D h_tof1("h_tof1", "ToF Eta1", 12, -6, 6);
    TH1D h_tof2("h_tof2", "ToF Eta2", 12, -6, 6);

    int event_count = 0;
    while (reader_track.Next()) {
        event_count++;
        int n_hits_eta1 = tof1->size();
        std::cout << "Event " << event_count << ": ToF Eta1 hits = " << n_hits_eta1 << std::endl;
        if (tof1.GetSetupStatus() == 0) {
            for (int t : *tof1) {
                h_tof1.Fill(t);
                std::cout << "ToF Eta1: " << t << std::endl;
                tof_results.time_of_flight_eta1.push_back(t);
            }
        }
        if (tof2.GetSetupStatus() == 0) {
            for (int t : *tof2) {
                h_tof2.Fill(t);
                std::cout << "ToF Eta2: " << t << std::endl;
                tof_results.time_of_flight_eta2.push_back(t);
            }
        }
    }

    // Reusable Lambda to perform the Two-Pass Gaussian Fit
    auto fitAndExtract = [](TH1D& hist, double& mean_out, ErrorRange& mean_err_out,
                            double& res_out, ErrorRange& res_err_out) {

        if (hist.GetEntries() < 20) return;

        // Pass 1: Global fit to find the general location of the peak
        TFitResultPtr r1 = hist.Fit("gaus", "Q0S");

        if (static_cast<int>(r1) == 0) { // If Fit 1 succeeded
            double p_mean  = r1->Parameter(1);
            double p_sigma = r1->Parameter(2);

            // Pass 2: Core fit restricted to +/- 1.5 sigma to ignore non-Gaussian tails
            TFitResultPtr r2 = hist.Fit("gaus", "Q0S", "", p_mean - 1.5 * p_sigma, p_mean + 1.5 * p_sigma);

            if (static_cast<int>(r2) == 0) {
                mean_out = r2->Parameter(1);
                mean_err_out = ErrorRange{r2->Error(1)};

                res_out = r2->Parameter(2) / std::sqrt(2.0);
                res_err_out = ErrorRange{r2->Error(2) / std::sqrt(2.0)};
            }
        }
    };

    fitAndExtract(h_tof1,
                  tof_results.avg_time_of_flight_eta1, tof_results.avg_time_of_flight_eta1_error,
                  time_resolution_results.time_resolution_eta1, time_resolution_results.time_resolution_eta1_error);

    fitAndExtract(h_tof2,
                  tof_results.avg_time_of_flight_eta2, tof_results.avg_time_of_flight_eta2_error,
                  time_resolution_results.time_resolution_eta2, time_resolution_results.time_resolution_eta2_error);
}

}

// ==========================================================================================
// DataAnalyzer class implementation for analyzing processed DCT data and plotting results
// ==========================================================================================
DataAnalyzer::DataAnalyzer(const std::string& config_file_path, const std::string& output_directory_path)
    : _config_path({config_file_path})
    , _output_directory(output_directory_path) {
}

DataAnalyzer::~DataAnalyzer() {
}

void DataAnalyzer::producePerFileStats(TFile* input_file) {

    // Set up output directories for per-file statistics and plots
    std::vector<std::string> dir_names = {
        "analysis/strip",
        "analysis/tot",
        "analysis/dt_strip",
        "analysis/tot_strip",
        "analysis/multiplicity_strip",
        "analysis/delay_strip",
        "analysis/tof"
    };
    PathUtils::setupDirectories(input_file, dir_names);

    // Produce relevant plots for this file using helper functions
    perFileHelpers::plotStrip(input_file);
    perFileHelpers::plotToT(input_file);
    perFileHelpers::plotDtVsStrip(input_file);
    perFileHelpers::plotToTVsStrip(input_file);
    perFileHelpers::plotMultiplicityAndDelayVsStrip(input_file);
    perFileHelpers::plotAdjacentToFs(input_file);

    // Write analysis directory into the file
    input_file->Write("", TObject::kOverwrite);

    // Export per-file plots to PDF
    std::string root_file_path = input_file->GetName();
    std::string target_plots_dir = (_output_directory / "plots" / std::filesystem::path(input_file->GetName()).stem()).string();
    PlotterHelpers::BatchExporter::autoExportToATLASPDF(root_file_path, target_plots_dir);
}

void DataAnalyzer::produceSummaryStats() {

    using namespace summaryHelpers;

    // Process each config file and build the list of measurement entries and summaries
    std::cout << _output_directory << std::endl;
    std::filesystem::create_directories(_output_directory / "root_summaries");

    // Prepare a measurement structure (measurement's metadata + data/statistics) and a root
    // file for this config
    ScanData scan;
    scan.config_path = _config_path;
    scan.metadata = ConfigUtils::parseMeasurementMetadata(_config_path);

    std::filesystem::path config_stem = std::filesystem::path(_config_path).stem();
    std::filesystem::path summary_root_path = _output_directory / "root_summaries" / (config_stem.string() + "_summary.root");

    // Create a ROOT file and tree to store summary statistics for this config
    TFile summary_root_file(summary_root_path.string().c_str(), "RECREATE");
    TTree* summary_tree = new TTree("summary", "summary");

    MeasurementMetadata metadata;
    MeasurementData data;

    // Set up branches for the summary tree
    setupBranches(summary_tree, metadata, data);

    // Process each measurement entry in this config file
    for (const auto& metadata_entry : scan.metadata) {
        if (metadata_entry.root_file.empty()) {
            std::cout << "Warning: Skipping entry '" << metadata_entry.name << "' due to missing ROOT file path." << std::endl;
            continue;
        }

        TFile* input_file = TFile::Open(metadata_entry.root_file.c_str(), "UPDATE");
        if (!input_file || input_file->IsZombie()) {
            std::cout << "Error: Failed to open ROOT file '" << metadata_entry.root_file << "' for entry '" << metadata_entry.name << "'. Skipping this entry." << std::endl;
            if (input_file) {
                std::cout << "Closing and deleting invalid ROOT file object for entry '" << metadata_entry.name << "'." << std::endl;
                input_file->Close();
                delete input_file;
            }
            continue;
        }

        data.clear();               // Reset data for this entry
        metadata = metadata_entry;  // Copy metadata for this entry

        // Calculate and fill per-file relevant statistics for this measurement entry and save into the input ROOT file
        producePerFileStats(input_file);

        // ----------------------------------------------------------------------------------
        // Extract efficiency values from histograms for this measurement entry
        for (int layer = 0; layer < LAYER_COUNT; ++layer) {

            // Helper lambda to safely read points and asymmetric errors out of a TGraphAsymmErrors object
            auto extractGraphData = [](TGraphAsymmErrors* graph, double effs[4], ErrorRange errs[4]) {
                for (int i = 0; i < 4; ++i) {
                    double dummy_x = 0.0;
                    // GetPoint stores the Y value directly inside the target efficiency index array
                    graph->GetPoint(i, dummy_x, effs[i]);
                    
                    // Extract the separate low and high bounds for the ErrorRange struct
                    errs[i].low  = graph->GetErrorYlow(i);
                    errs[i].high = graph->GetErrorYhigh(i);
                }
            };

            // ------------------------------------------------------------------------------
            // 1. External trigger only
            std::string path_ext = "efficiency_graphs/external_trigger/eff_external_trigger_layer" + std::to_string(layer);
            TGraphAsymmErrors* graph_ext = dynamic_cast<TGraphAsymmErrors*>(input_file->Get(path_ext.c_str()));
            if (graph_ext) {
                double effs[4];
                ErrorRange errs[4];
                extractGraphData(graph_ext, effs, errs);

                data.efficiency_results.eta1_efficiency_external[layer] = effs[0];
                data.efficiency_results.eta2_efficiency_external[layer] = effs[1];
                data.efficiency_results.eta_or_efficiency_external[layer] = effs[2];
                data.efficiency_results.eta_and_efficiency_external[layer] = effs[3];

                data.efficiency_results.eta1_efficiency_external_error[layer] = errs[0];
                data.efficiency_results.eta2_efficiency_external_error[layer] = errs[1];
                data.efficiency_results.eta_or_efficiency_external_error[layer] = errs[2];
                data.efficiency_results.eta_and_efficiency_external_error[layer] = errs[3];
            }

            // ------------------------------------------------------------------------------
            // 2. External trigger + RPC
            std::string path_rpc = "efficiency_graphs/external_plus_rpc_trigger/eff_rpc_layer" + std::to_string(layer);
            TGraphAsymmErrors* graph_rpc = dynamic_cast<TGraphAsymmErrors*>(input_file->Get(path_rpc.c_str()));
            if (graph_rpc) {
                double effs[4];
                ErrorRange errs[4];
                extractGraphData(graph_rpc, effs, errs);

                data.efficiency_results.eta1_efficiency_rpc[layer] = effs[0];
                data.efficiency_results.eta2_efficiency_rpc[layer] = effs[1];
                data.efficiency_results.eta_or_efficiency_rpc[layer] = effs[2];
                data.efficiency_results.eta_and_efficiency_rpc[layer] = effs[3];

                data.efficiency_results.eta1_efficiency_rpc_error[layer] = errs[0];
                data.efficiency_results.eta2_efficiency_rpc_error[layer] = errs[1];
                data.efficiency_results.eta_or_efficiency_rpc_error[layer] = errs[2];
                data.efficiency_results.eta_and_efficiency_rpc_error[layer] = errs[3];
            }

            // ------------------------------------------------------------------------------
            // 3. Track external trigger only
            std::string path_track = "efficiency_graphs/track_external_trigger/track_eff_external_trigger_layer" + std::to_string(layer);
            TGraphAsymmErrors* graph_track = dynamic_cast<TGraphAsymmErrors*>(input_file->Get(path_track.c_str()));
            if (graph_track) {
                double effs[4];
                ErrorRange errs[4];
                extractGraphData(graph_track, effs, errs);

                data.efficiency_results_tracks.eta1_efficiency_external[layer] = effs[0];
                data.efficiency_results_tracks.eta2_efficiency_external[layer] = effs[1];
                data.efficiency_results_tracks.eta_or_efficiency_external[layer] = effs[2];
                data.efficiency_results_tracks.eta_and_efficiency_external[layer] = effs[3];

                data.efficiency_results_tracks.eta1_efficiency_external_error[layer] = errs[0];
                data.efficiency_results_tracks.eta2_efficiency_external_error[layer] = errs[1];
                data.efficiency_results_tracks.eta_or_efficiency_external_error[layer] = errs[2];
                data.efficiency_results_tracks.eta_and_efficiency_external_error[layer] = errs[3];
            }

            // ------------------------------------------------------------------------------
            // 4. Track external trigger + RPC
            std::string path_track_rpc = "efficiency_graphs/track_external_plus_rpc_trigger/track_eff_rpc_layer" + std::to_string(layer);
            TGraphAsymmErrors* graph_track_rpc = dynamic_cast<TGraphAsymmErrors*>(input_file->Get(path_track_rpc.c_str()));
            if (graph_track_rpc) {
                double effs[4];
                ErrorRange errs[4];
                extractGraphData(graph_track_rpc, effs, errs);

                data.efficiency_results_tracks.eta1_efficiency_rpc[layer] = effs[0];
                data.efficiency_results_tracks.eta2_efficiency_rpc[layer] = effs[1];
                data.efficiency_results_tracks.eta_or_efficiency_rpc[layer] = effs[2];
                data.efficiency_results_tracks.eta_and_efficiency_rpc[layer] = effs[3];

                data.efficiency_results_tracks.eta1_efficiency_rpc_error[layer] = errs[0];
                data.efficiency_results_tracks.eta2_efficiency_rpc_error[layer] = errs[1];
                data.efficiency_results_tracks.eta_or_efficiency_rpc_error[layer] = errs[2];
                data.efficiency_results_tracks.eta_and_efficiency_rpc_error[layer] = errs[3];
            }
        }

        // ----------------------------------------------------------------------------------
        {  // Average cluster size calculation
        TTree* cluster_tree = dynamic_cast<TTree*>(input_file->Get("Clusterization"));
        if (!cluster_tree) {
            std::cerr << "Error: Clusterization tree not found in file " << input_file->GetName() << std::endl;
            std::exit(EXIT_FAILURE);
        }
        TTreeReader reader(cluster_tree);
        TTreeReaderValue<std::vector<int>> cluster_eta1(reader, "cluster_size_eta1");
        TTreeReaderValue<std::vector<int>> cluster_eta2(reader, "cluster_size_eta2");
        TTreeReaderValue<std::vector<int>> cluster_eta1_layer0(reader, "cluster_size_eta1_layer0");
        TTreeReaderValue<std::vector<int>> cluster_eta1_layer1(reader, "cluster_size_eta1_layer1");
        TTreeReaderValue<std::vector<int>> cluster_eta1_layer2(reader, "cluster_size_eta1_layer2");
        TTreeReaderValue<std::vector<int>> cluster_eta2_layer0(reader, "cluster_size_eta2_layer0");
        TTreeReaderValue<std::vector<int>> cluster_eta2_layer1(reader, "cluster_size_eta2_layer1");
        TTreeReaderValue<std::vector<int>> cluster_eta2_layer2(reader, "cluster_size_eta2_layer2");

        long long total_clusters_eta1 = 0;
        long long total_clusters_eta2 = 0;
        long long total_cluster_size_eta1 = 0;
        long long total_cluster_size_eta2 = 0;

        long long sum_sq_cluster_size_eta1 = 0;
        long long sum_sq_cluster_size_eta2 = 0;

        long long total_clusters_eta1_layers[3] = {0, 0, 0};
        long long total_clusters_eta2_layers[3] = {0, 0, 0};
        long long total_cluster_size_eta1_layers[3] = {0, 0, 0};
        long long total_cluster_size_eta2_layers[3] = {0, 0, 0};

        long long sum_sq_cluster_size_eta1_layers[3] = {0, 0, 0};
        long long sum_sq_cluster_size_eta2_layers[3] = {0, 0, 0};

        while (reader.Next()) {
            // Define a generic, inline accumulator helper
            auto accumulateBranch = [](auto& reader_val, long long& count, long long& sum, long long& sum_sq) {
                if (reader_val.GetSetupStatus() == 0) {
                    count += static_cast<long long>(reader_val->size());
                    for (int size : *reader_val) {
                        sum += size;
                        sum_sq += (static_cast<long long>(size) * size);
                    }
                }
            };

            // Accumulate global metrics
            accumulateBranch(cluster_eta1, total_clusters_eta1, total_cluster_size_eta1, sum_sq_cluster_size_eta1);
            accumulateBranch(cluster_eta2, total_clusters_eta2, total_cluster_size_eta2, sum_sq_cluster_size_eta2);

            // Accumulate layer-by-layer metrics
            accumulateBranch(cluster_eta1_layer0, total_clusters_eta1_layers[0], total_cluster_size_eta1_layers[0], sum_sq_cluster_size_eta1_layers[0]);
            accumulateBranch(cluster_eta1_layer1, total_clusters_eta1_layers[1], total_cluster_size_eta1_layers[1], sum_sq_cluster_size_eta1_layers[1]);
            accumulateBranch(cluster_eta1_layer2, total_clusters_eta1_layers[2], total_cluster_size_eta1_layers[2], sum_sq_cluster_size_eta1_layers[2]);
            
            accumulateBranch(cluster_eta2_layer0, total_clusters_eta2_layers[0], total_cluster_size_eta2_layers[0], sum_sq_cluster_size_eta2_layers[0]);
            accumulateBranch(cluster_eta2_layer1, total_clusters_eta2_layers[1], total_cluster_size_eta2_layers[1], sum_sq_cluster_size_eta2_layers[1]);
            accumulateBranch(cluster_eta2_layer2, total_clusters_eta2_layers[2], total_cluster_size_eta2_layers[2], sum_sq_cluster_size_eta2_layers[2]);
        }

        // Helper lambda to calculate mean and standard error of the mean (SEM)
        auto calcMeanAndError = [](long long N, long long sum, long long sum_sq, double& mean_out, ErrorRange& err_out) {
            if (N <= 0) {
                mean_out = 0.0;
                err_out = ErrorRange{0.0};
                return;
            }
            mean_out = static_cast<double>(sum) / N;
            
            if (N > 1) {
                double variance = (static_cast<double>(sum_sq) - (static_cast<double>(sum) * sum) / N) / (N - 1);
                if (variance < 0.0) variance = 0.0; 
                double std_dev = std::sqrt(variance);
                err_out = ErrorRange{std_dev / std::sqrt(N)};
            } else {
                err_out = ErrorRange{0.0}; // Cannot calculate error with only one cluster
            }
        };

        // Calculate global metrics
        calcMeanAndError(total_clusters_eta1, total_cluster_size_eta1, sum_sq_cluster_size_eta1, 
                         data.cluster_size_results.avg_cluster_size_eta1, 
                         data.cluster_size_results.avg_cluster_size_eta1_error);

        calcMeanAndError(total_clusters_eta2, total_cluster_size_eta2, sum_sq_cluster_size_eta2, 
                         data.cluster_size_results.avg_cluster_size_eta2, 
                         data.cluster_size_results.avg_cluster_size_eta2_error);

        // Calculate layer-by-layer metrics
        for (int layer_idx = 0; layer_idx < 3; ++layer_idx) {
            calcMeanAndError(total_clusters_eta1_layers[layer_idx], 
                             total_cluster_size_eta1_layers[layer_idx], 
                             sum_sq_cluster_size_eta1_layers[layer_idx],
                             data.cluster_size_results.avg_cluster_size_eta1_layers[layer_idx],
                             data.cluster_size_results.avg_cluster_size_eta1_layers_error[layer_idx]);

            calcMeanAndError(total_clusters_eta2_layers[layer_idx], 
                             total_cluster_size_eta2_layers[layer_idx], 
                             sum_sq_cluster_size_eta2_layers[layer_idx],
                             data.cluster_size_results.avg_cluster_size_eta2_layers[layer_idx],
                             data.cluster_size_results.avg_cluster_size_eta2_layers_error[layer_idx]);
        }
        }   // Close the scope for cluster size calculation to avoid variable name conflicts

        // ----------------------------------------------------------------------------------
        // Noise rate calculation
        summaryHelpers::getRate(input_file, data.noise_rate_results);

        // ----------------------------------------------------------------------------------
        // Average Time over Threshold (ToT) calculation
        summaryHelpers::getAverageToT(input_file, data.tot_results, false);
        summaryHelpers::getAverageToT(input_file, data.tot_results_tracks, true);

        // ----------------------------------------------------------------------------------
        // Average multiplicity calculation
        summaryHelpers::getAverageMultiplicity(input_file, data.multiplicity_results, false);
        summaryHelpers::getAverageMultiplicity(input_file, data.multiplicity_results_tracks, true);

        // ----------------------------------------------------------------------------------
        // Average Time of Flight (ToF) calculation
        summaryHelpers::processToF(input_file, data.tof_results, data.time_resolution_results);

        // Fill the summary tree with the extracted statistics for this measurement entry
        summary_tree->Fill();
        scan.data.push_back(data);

        // Clean up and close the input ROOT file for this measurement entry
        input_file->Close();
        delete input_file;
    }

    summary_root_file.cd();
    summary_tree->Write("", TObject::kOverwrite);
    summary_root_file.Close();
}
