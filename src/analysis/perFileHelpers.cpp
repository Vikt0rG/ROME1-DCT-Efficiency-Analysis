#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TProfile2D.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TLegend.h>

#include "utils.hpp"
#include "constants.hpp"

#include "core/hit.hpp"
#include "core/event.hpp"
#include "core/cluster.hpp"
#include "core/track.hpp"
#include "analysis/perFileHelpers.hpp"


namespace perFileHelpers {

void plotEfficiencyMap(TFile* input_file, bool use_external_trigger = false) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* eff_dir = analysis_dir->GetDirectory("efficiency_map");
    eff_dir->cd();

    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));

    requireEqualTreeEntries({
        {"ProcessedData", proc_tree},
        {"TrackReconstruction", track_tree}
    });

    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);

    TTreeReaderValue<int> proc_trigger_time(readerProcData, "proc_trigger_time");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> dts(readerProcData, "proc_dt_time1_time2");
    TTreeReaderValue<std::vector<int>> proc_dt_time1_trigger(readerProcData, "proc_dt_time1_trigger");
    TTreeReaderValue<std::vector<int>> proc_dt_time2_trigger(readerProcData, "proc_dt_time2_trigger");
    TTreeReaderValue<std::vector<bool>> has_valid_time_eta1(readerProcData, "has_valid_time_eta1");
    TTreeReaderValue<std::vector<bool>> has_valid_time_eta2(readerProcData, "has_valid_time_eta2");

    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    const int n_bins_x = STRIPS_PER_LAYER;

    const int n_bins_y = N_BINS_DT;
    const int dt_max_ticks = n_bins_y / 2;
    const int dt_min_ticks = -dt_max_ticks;
    const double dt_range_ticks = dt_max_ticks - dt_min_ticks;

    // Standard raw array allocation matching dynamic n_bins_y
    std::vector<std::vector<std::vector<int>>> passed_counts(LAYER_COUNT, std::vector<std::vector<int>>(STRIPS_PER_LAYER, std::vector<int>(n_bins_y, 0)));
    std::vector<std::vector<std::vector<int>>> total_counts(LAYER_COUNT, std::vector<std::vector<int>>(STRIPS_PER_LAYER, std::vector<int>(n_bins_y, 0)));

    std::vector<std::vector<std::vector<int>>> track_passed_counts(LAYER_COUNT, std::vector<std::vector<int>>(STRIPS_PER_LAYER, std::vector<int>(n_bins_y, 0)));
    std::vector<std::vector<std::vector<int>>> track_total_counts(LAYER_COUNT, std::vector<std::vector<int>>(STRIPS_PER_LAYER, std::vector<int>(n_bins_y, 0)));

    int event_idx = 0;

    while (readerProcData.Next() && readerTrackData.Next()) {
        event_idx++;
        if (use_external_trigger && *proc_trigger_time == -1) continue;

        try {
            requireEqualDataSizes({
                {"proc_strip", strips->size()},
                {"proc_layer", layers->size()},
                {"proc_dt_time1_time2", dts->size()},
                {"proc_dt_time1_trigger", proc_dt_time1_trigger->size()},
                {"proc_dt_time2_trigger", proc_dt_time2_trigger->size()},
                {"has_valid_time_eta1", has_valid_time_eta1->size()},
                {"has_valid_time_eta2", has_valid_time_eta2->size()},
                {"in_valid_track_eta1", in_valid_track_eta1->size()},
                {"in_valid_track_eta2", in_valid_track_eta2->size()}
            });
        } catch (const std::runtime_error& e) {
            std::cerr << "Error in event #" << event_idx << ": " << e.what() << std::endl;
            throw;
        }

        std::vector<std::vector<std::vector<bool>>> cell_hits(LAYER_COUNT, std::vector<std::vector<bool>>(STRIPS_PER_LAYER, std::vector<bool>(n_bins_y, false)));
        std::vector<std::vector<std::vector<bool>>> track_cell_hits(LAYER_COUNT, std::vector<std::vector<bool>>(STRIPS_PER_LAYER, std::vector<bool>(n_bins_y, false)));

        const size_t n_hits = strips->size();

        for (size_t i = 0; i < n_hits; ++i) {

            if (!(*has_valid_time_eta1)[i] || !(*has_valid_time_eta2)[i]) continue;

            int layer = (*layers)[i];
            if (layer < 0 || layer >= LAYER_COUNT) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            if (strip < 0 || strip >= STRIPS_PER_LAYER) continue;

            int dt1_trig = (*proc_dt_time1_trigger)[i];
            int dt2_trig = (*proc_dt_time2_trigger)[i];

            // Time-window cut w.r.t. trigger time
            if (dt1_trig >= DEFAULT_DT_MIN && dt1_trig <= DEFAULT_DT_MAX &&
                dt2_trig >= DEFAULT_DT_MIN && dt2_trig <= DEFAULT_DT_MAX) 
            {
                double dt_ticks = static_cast<double>((*dts)[i]);

                // Linear Mapping Formula: Maps signed values to [0, n_bins_y - 1]
                if (dt_ticks >= dt_min_ticks && dt_ticks < dt_max_ticks) {
                    int dt_idx = static_cast<int>(((dt_ticks - dt_min_ticks) / dt_range_ticks) * n_bins_y);
                    dt_idx = std::clamp(dt_idx, 0, n_bins_y - 1);

                    cell_hits[layer][strip][dt_idx] = true;

                    if ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]) {
                        track_cell_hits[layer][strip][dt_idx] = true;
                    }
                }
            }
        }

        // Tag-and-Probe Coincidence Computation across the coarse cells
        for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {
            for (int dt = 0; dt < n_bins_y; ++dt) {
                // Raw cell hits
                bool h0 = cell_hits[0][strip][dt];
                bool h1 = cell_hits[1][strip][dt];
                bool h2 = cell_hits[2][strip][dt];

                if (h1 && h2) {
                    total_counts[0][strip][dt]++;
                    if (h0) passed_counts[0][strip][dt]++;
                }
                if (h0 && h2) {
                    total_counts[1][strip][dt]++;
                    if (h1) passed_counts[1][strip][dt]++;
                }
                if (h0 && h1) {
                    total_counts[2][strip][dt]++;
                    if (h2) passed_counts[2][strip][dt]++;
                }

                // Track cell hits
                bool th0 = track_cell_hits[0][strip][dt];
                bool th1 = track_cell_hits[1][strip][dt];
                bool th2 = track_cell_hits[2][strip][dt];

                if (th1 && th2) {
                    track_total_counts[0][strip][dt]++;
                    if (th0) track_passed_counts[0][strip][dt]++;
                }
                if (th0 && th2) {
                    track_total_counts[1][strip][dt]++;
                    if (th1) track_passed_counts[1][strip][dt]++;
                }
                if (th0 && th1) {
                    track_total_counts[2][strip][dt]++;
                    if (th2) track_passed_counts[2][strip][dt]++;
                }
            }
        }
    }

    // Physical bounds in cm centered at 0
    const double x_min_cm = 0.0;
    const double x_max_cm = STRIPS_PER_LAYER * STRIP_WIDTH_CM;

    const double y_min_cm = (SIGNAL_SPEED * dt_min_ticks * TIME_TICK_NS) / 2.0;
    const double y_max_cm = (SIGNAL_SPEED * dt_max_ticks * TIME_TICK_NS) / 2.0;

    struct MapConfig {
        std::string suffix;
        std::string title_tag;
        const std::vector<std::vector<std::vector<int>>>& passed_ptr;
        const std::vector<std::vector<std::vector<int>>>& total_ptr;
    };

    std::vector<MapConfig> configs = {
        {"before_reco", "Before Track Reco", passed_counts, total_counts},
        {"after_reco", "After Track Reco", track_passed_counts, track_total_counts}
    };

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    for (const auto& cfg : configs) {
        for (int layer = 0; layer < LAYER_COUNT; ++layer) {
            TH2F* h_eff_map = new TH2F(
                Form("h2d_efficiency_map_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    cfg.suffix.c_str(),
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv),
                Form("Efficiency Map: Layer %d: %s;X [cm];Y [cm];Efficiency", layer, cfg.title_tag.c_str()),
                n_bins_x, x_min_cm, x_max_cm,
                n_bins_y, y_min_cm, y_max_cm
            );

            h_eff_map->SetMinimum(0.0);
            h_eff_map->SetMaximum(1.0);

            for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {
                for (int dt = 0; dt < n_bins_y; ++dt) {
                    int passed = cfg.passed_ptr[layer][strip][dt];
                    int total = cfg.total_ptr[layer][strip][dt];

                    int bx = strip + 1;
                    int by = dt + 1;

                    if (total > 0) {
                        double eff = static_cast<double>(passed) / static_cast<double>(total);
                        double err = std::sqrt(eff * (1.0 - eff) / static_cast<double>(total));

                        h_eff_map->SetBinContent(bx, by, eff);
                        h_eff_map->SetBinError(bx, by, err);
                    } else {
                        h_eff_map->SetBinContent(bx, by, 0.0);
                        h_eff_map->SetBinError(bx, by, 0.0);
                    }
                }
            }

            h_eff_map->Write("", TObject::kOverwrite);
            delete h_eff_map;
        }
    }
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
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid track reconstruction tree for analysis." << std::endl;
        return;
    }

    // Read processed data into vectors and create strip distributions for each layer and side
    TTreeReader readerInputData(input_data_tree);
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> raw_time1(readerInputData, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>> raw_time2(readerInputData, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    const int nConfigs = 4;
    const char* categories[nConfigs] = {
        "strip_eta1_after_reco", "strip_eta2_after_reco",
        "strip_eta1_rejected", "strip_eta2_rejected"
    };

    std::map<std::string, std::map<int, TH1*>> strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH1F(Form("h1d_%s_layer%d", categories[c], layer),
                                  Form("Layer %d;Strip;Hits", layer),
                                  STRIPS_PER_LAYER, 0, STRIPS_PER_LAYER);
            strip_histograms[categories[c]][layer] = hist;
        }
    }

    while (readerInputData.Next() && readerProcData.Next() && readerTrackData.Next()) {
        const size_t n_hits = std::min({strips->size(), layers->size(), 
                                        raw_time1->size(), raw_time2->size(), 
                                        in_valid_track_eta1->size(), in_valid_track_eta2->size()});

        for (size_t i = 0; i < n_hits; ++i) {
            int layer = (*layers)[i];
            if (layer < 0 || layer > 2) continue;

            int strip = remapStrip((*strips)[i]);

            if ((*raw_time1)[i] != 0) {
                if ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]) {
                    strip_histograms["strip_eta1_after_reco"][layer]->Fill(strip);
                } else {
                    strip_histograms["strip_eta1_rejected"][layer]->Fill(strip);
                }
            }
            if ((*raw_time2)[i] != 0) {
                if ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]) {
                    strip_histograms["strip_eta2_after_reco"][layer]->Fill(strip);
                } else {
                    strip_histograms["strip_eta2_rejected"][layer]->Fill(strip);
                }
            }
        }
    }

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            strip_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
        }
    }

    struct StackPairing {
        std::string side;
        std::pair<std::string, std::string> pair;
        std::string suffix;
    };

    std::vector<StackPairing> stack_pairings = {
        {"eta1", {"strip_eta1_after_reco", "strip_eta1_rejected"}, "reco_stack"},
        {"eta2", {"strip_eta2_after_reco", "strip_eta2_rejected"}, "reco_stack"}
    };

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    for (const auto& pairing : stack_pairings) {
        for (int layer : {0, 1, 2}) {
            THStack* stack = new THStack(
                Form("h1d_strip_%s_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    pairing.side.c_str(),
                    pairing.suffix.c_str(),
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv),
                Form("Layer %d;Strip;Hits", layer)
            );

            TH1* after_hist = strip_histograms[pairing.pair.first][layer];
            TH1* rejected_hist = strip_histograms[pairing.pair.second][layer];

            after_hist->SetFillColor(kRed + 1);
            rejected_hist->SetFillColor(kAzure - 3);

            stack->Add(after_hist);
            stack->Add(rejected_hist);

            stack->Write("", TObject::kOverwrite);
            delete stack;
        }
    }

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            delete strip_histograms[categories[c]][layer];
        }
    }
}

void extractBeamSpotRoI(TH2* agg_hist, RoI& roi, double core_frac, double halo_frac) {
    roi.clear();
    if (!agg_hist) return;

    double max_hits = agg_hist->GetMaximum();
    if (max_hits <= 0) return;

    double core_thresh = max_hits * core_frac;
    double halo_thresh = max_hits * halo_frac;

    for (int bx = 1; bx <= agg_hist->GetNbinsX(); ++bx) {
        for (int by = 1; by <= agg_hist->GetNbinsY(); ++by) {

            double content = agg_hist->GetBinContent(bx, by);
            int strip = std::round(agg_hist->GetXaxis()->GetBinCenter(bx));
            int dt = std::round(agg_hist->GetYaxis()->GetBinCenter(by));

            if (content > core_thresh) {
                roi.beam_region.push_back({strip, dt});
            } else if (content > halo_thresh) {
                roi.halo_region.push_back({strip, dt});
            }
        }
    }
}

void plotDtVsStrip(TFile* input_file, RoI& region_of_interest) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* dt_strip_dir = analysis_dir->GetDirectory("dt_strip");
    dt_strip_dir->cd();

    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!proc_tree || proc_tree->IsZombie() || !track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid processed or track tree for analysis." << std::endl;
        return;
    }

    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> dts(readerProcData, "proc_dt_time1_time2");
    TTreeReaderValue<std::vector<bool>> has_valid_time_eta1(readerProcData, "has_valid_time_eta1");
    TTreeReaderValue<std::vector<bool>> has_valid_time_eta2(readerProcData, "has_valid_time_eta2");

    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    const int nConfigs = 2;
    const char* reco_status[nConfigs] = {"before_reco", "after_reco"};
    const char* comments[nConfigs] = {"Before Track Reco", "After Track Reco"};

    const int n_bins_x = STRIPS_PER_LAYER;
    const double x_min_strip = 0.0;
    const double x_max_strip = static_cast<double>(STRIPS_PER_LAYER);

    const int n_bins_y = N_BINS_DT;
    const int dt_max_ticks = n_bins_y / 2;
    const int dt_min_ticks = -dt_max_ticks;

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    std::map<std::string, std::map<int, TH2F*>> dt_strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(
                Form("h2d_dt_strip_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    reco_status[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Strip;#Delta#it{t} [Ticks];Entries", layer, comments[c]),
                n_bins_x, x_min_strip, x_max_strip,
                n_bins_y, dt_min_ticks, dt_max_ticks
            );



            dt_strip_histograms[reco_status[c]][layer] = hist;
        }
    }

    long long n_events = 0;
    while (readerProcData.Next() && readerTrackData.Next()) {
        n_events++;

        const size_t n_hits = std::min({strips->size(), layers->size(), dts->size()});
        for (size_t i = 0; i < n_hits; ++i) {
            if (!(*has_valid_time_eta1)[i] || !(*has_valid_time_eta2)[i]) continue;

            int layer = (*layers)[i];
            if (layer < 0 || layer >= LAYER_COUNT) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            if (strip < 0 || strip >= STRIPS_PER_LAYER) continue;

            double dt_ticks = static_cast<double>((*dts)[i]);

            dt_strip_histograms["before_reco"][layer]->Fill(strip, dt_ticks);
            if ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]) {
                dt_strip_histograms["after_reco"][layer]->Fill(strip, dt_ticks);
            }
        }
    }

    TH2F* h_agg = (TH2F*)dt_strip_histograms["after_reco"][0]->Clone("h2d_agg_tracks");
    h_agg->Add(dt_strip_histograms["after_reco"][1]);
    h_agg->Add(dt_strip_histograms["after_reco"][2]);
    h_agg->SetTitle("Aggregated dt vs Strip for All Layers;Strip;#Delta#it{t} [Ticks]");

    double beam_threshold = 0.50;
    double halo_threshold = 0.05;
    extractBeamSpotRoI(h_agg, region_of_interest, beam_threshold, halo_threshold);

    TH2F* h_roi_map = new TH2F(
        "h2d_roi_map", "Beam Spot Region of Interest;Strip;#Delta#it{t} [Ticks]",
        n_bins_x, x_min_strip, x_max_strip, n_bins_y, dt_min_ticks, dt_max_ticks
    );

    // Default 1.0 (No Beam)
    for (int bx = 1; bx <= h_roi_map->GetNbinsX(); ++bx) {
        for (int by = 1; by <= h_roi_map->GetNbinsY(); ++by) {
            h_roi_map->SetBinContent(bx, by, 1.0);
        }
    }

    // Mark Halo region (2.0)
    for (const auto& pt : region_of_interest.halo_region) {
        int bx = h_roi_map->GetXaxis()->FindBin(pt.first + 0.5);
        int by = h_roi_map->GetYaxis()->FindBin(pt.second);
        h_roi_map->SetBinContent(bx, by, 2.0);
    }

    // Mark Beam region (3.0)
    for (const auto& pt : region_of_interest.beam_region) {
        int bx = h_roi_map->GetXaxis()->FindBin(pt.first + 0.5);
        int by = h_roi_map->GetYaxis()->FindBin(pt.second);
        h_roi_map->SetBinContent(bx, by, 3.0);
    }

    h_agg->Write("", TObject::kOverwrite);
    h_roi_map->Write("", TObject::kOverwrite);

    delete h_agg;
    delete h_roi_map;

    double dx_cm = STRIP_WIDTH_CM;
    double dt_bin_width_ns = TIME_TICK_NS;
    double dy_cm = (SIGNAL_SPEED * dt_bin_width_ns) / 2.0;

    double cell_area_cm2 = dx_cm * dy_cm;
    double total_exposure = n_events * TRIGGER_TIME_WINDOW * cell_area_cm2;

    const double x_min_cm = 0.0;
    const double x_max_cm = STRIPS_PER_LAYER * STRIP_WIDTH_CM;
    const double y_min_cm = (SIGNAL_SPEED * dt_min_ticks * TIME_TICK_NS) / 2.0;
    const double y_max_cm = (SIGNAL_SPEED * dt_max_ticks * TIME_TICK_NS) / 2.0;

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            TH2F* h_count = dt_strip_histograms[reco_status[c]][layer];

            h_count->Write("", TObject::kOverwrite);

            if (total_exposure > 0.0) {

                std::string filter_str = Form("%.1f", _metadata.filter);
                std::replace(filter_str.begin(), filter_str.end(), '.', '_');
                TH2F* h_rate = new TH2F(
                    Form("h2d_rate_map_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                        reco_status[c],
                        layer,
                        _metadata.mixture.c_str(),
                        _metadata.lv_setting,
                        filter_str.c_str(),
                        _metadata.scanned_hv
                    ),
                    Form("Rate Map: Layer %d: %s;X [cm];Y [cm];Rate [Hz/cm^{2}]", layer, comments[c]),
                    n_bins_x, x_min_cm, x_max_cm,
                    n_bins_y, y_min_cm, y_max_cm
                );

                for (int bx = 1; bx <= n_bins_x; ++bx) {
                    for (int by = 1; by <= n_bins_y; ++by) {
                        double counts = h_count->GetBinContent(bx, by);
                        if (counts > 0.0) {
                            h_rate->SetBinContent(bx, by, counts / total_exposure);
                            h_rate->SetBinError(bx, by, std::sqrt(counts) / total_exposure);
                        } else {
                            h_rate->SetBinContent(bx, by, 0.0);
                            h_rate->SetBinError(bx, by, 0.0);
                        }
                    }
                }

                h_rate->Write("", TObject::kOverwrite);
                delete h_rate;
            }

            delete h_count;
        }
    }
}

void plotCS(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* cs_dir = analysis_dir->GetDirectory("cluster_size");
    if (!analysis_dir || !cs_dir) {
        std::cerr << "Error: Missing analysis or cluster size directory in input file." << std::endl;
        return;
    }
    cs_dir->cd();

    TTree* cluster_tree = dynamic_cast<TTree*>(input_file->Get("Clusterization"));
    if (!cluster_tree || cluster_tree->IsZombie()) {
        std::cerr << "Error: Invalid clusterization tree for analysis." << std::endl;
        return;
    }

    int primary_hv = extractPrimaryHV(input_file->GetName());

    // Read processed data into vectors
    TTreeReader readerClusterData(cluster_tree);
    TTreeReaderValue<std::vector<int>> cluster_eta1_layer0(readerClusterData, "cluster_size_eta1_layer0");
    TTreeReaderValue<std::vector<int>> cluster_eta1_layer1(readerClusterData, "cluster_size_eta1_layer1");
    TTreeReaderValue<std::vector<int>> cluster_eta1_layer2(readerClusterData, "cluster_size_eta1_layer2");
    TTreeReaderValue<std::vector<int>> cluster_eta2_layer0(readerClusterData, "cluster_size_eta2_layer0");
    TTreeReaderValue<std::vector<int>> cluster_eta2_layer1(readerClusterData, "cluster_size_eta2_layer1");
    TTreeReaderValue<std::vector<int>> cluster_eta2_layer2(readerClusterData, "cluster_size_eta2_layer2");

    const int nConfigs = 2;
    const char* side_context[nConfigs] = {"eta1", "eta2"};
    const char* comments[nConfigs] = {"Side #eta+", "Side #eta-"};

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    std::map<std::string, std::map<int, TH1*>> strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            std::string hist_name = Form(
                "h1d_cs_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                side_context[c],
                layer,
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            );
            std::string hist_title = Form("Layer %d: %s (HV: %d V);Cluster Size [Strips];Counts",
                                          layer, comments[c], primary_hv);

            auto* hist = new TH1F(hist_name.c_str(), hist_title.c_str(), 6, 0.5, 6.5);

            strip_histograms[side_context[c]][layer] = hist;
        }
    }

    while (readerClusterData.Next()) {
        for (int cs : *cluster_eta1_layer0) strip_histograms["eta1"][0]->Fill(cs);
        for (int cs : *cluster_eta1_layer1) strip_histograms["eta1"][1]->Fill(cs);
        for (int cs : *cluster_eta1_layer2) strip_histograms["eta1"][2]->Fill(cs);

        for (int cs : *cluster_eta2_layer0) strip_histograms["eta2"][0]->Fill(cs);
        for (int cs : *cluster_eta2_layer1) strip_histograms["eta2"][1]->Fill(cs);
        for (int cs : *cluster_eta2_layer2) strip_histograms["eta2"][2]->Fill(cs);
    }

    cs_dir->cd();

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            TH1* hist = strip_histograms[side_context[c]][layer];

            if (hist->GetEntries() < 5) {
                hist->Write("", TObject::kOverwrite);
                continue;
            }

            // Shifted Poisson: x-1 represents the "extra" strips beyond the guaranteed 1.
            // [0] = Amplitude, [1] = Mean number of extra strips
            TF1* cs_fit = new TF1("cs_poisson", "[0] * TMath::Poisson(x - 1.0, [1])", 0.5, 6.5);

            double max_val = hist->GetMaximum();
            double raw_mean = hist->GetMean();

            cs_fit->SetParameter(0, max_val);

            double initial_lambda = std::max(raw_mean - 1.0, 0.01);
            cs_fit->SetParameter(1, initial_lambda);
            cs_fit->SetParLimits(1, 0.001, 5.0);

            cs_fit->SetLineColor(kRed);
            cs_fit->SetLineWidth(2);

            hist->Fit(cs_fit, "Q");

            hist->Write("", TObject::kOverwrite);
        }
    }
}

void plotToT(TFile* input_file, const RoI& region_of_interest) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* tot_dir = analysis_dir->GetDirectory("tot");
    tot_dir->cd();

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

    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> time1(readerProcData, "proc_time1");
    TTreeReaderValue<std::vector<int>> time2(readerProcData, "proc_time2");
    TTreeReaderValue<std::vector<int>> tot1(readerProcData, "proc_tot1");
    TTreeReaderValue<std::vector<int>> tot2(readerProcData, "proc_tot2");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> dts(readerProcData, "proc_dt_time1_time2");
    TTreeReaderValue<std::vector<bool>> has_valid_time_eta1(readerProcData, "has_valid_time_eta1");
    TTreeReaderValue<std::vector<bool>> has_valid_time_eta2(readerProcData, "has_valid_time_eta2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    const int nConfigs = 6;
    const char* categories[nConfigs] = {
        "tot_eta1_before_reco", "tot_eta2_before_reco",
        "tot_eta1_after_reco",  "tot_eta2_after_reco",
        "tot_eta1_beam",        "tot_eta2_beam"
    };
    const char* comments[nConfigs] = {
        "Before Track Reco", "Before Track Reco",
        "After Track Reco",  "After Track Reco",
        "Beam Spot Region",  "Beam Spot Region"
    };

    // 12 Multiplicity categories (primary, secondary hits per side)
    const int nMultiConfigs = 12;
    const char* multi_categories[nMultiConfigs] = {
        "tot_eta1_primary_before_reco", "tot_eta1_secondary_before_reco",
        "tot_eta2_primary_before_reco", "tot_eta2_secondary_before_reco",
        "tot_eta1_primary_after_reco",  "tot_eta1_secondary_after_reco",
        "tot_eta2_primary_after_reco",  "tot_eta2_secondary_after_reco",
        "tot_eta1_primary_beam",        "tot_eta1_secondary_beam",
        "tot_eta2_primary_beam",        "tot_eta2_secondary_beam"
    };

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    int max_expected_ticks = 150;
    const int nBins = max_expected_ticks;
    const double xMin = -0.5 * TIME_TICK_NS;
    const double xMax = (max_expected_ticks - 0.5) * TIME_TICK_NS;

    std::map<std::string, std::map<int, TH1*>> strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH1F(Form("h1d_%s_layer%d", categories[c], layer),
                            Form("Layer %d: %s;ToT [ns];Hits", layer, comments[c]),
                            nBins, xMin, xMax);
            strip_histograms[categories[c]][layer] = hist;
        }
    }

    for (int c = 0; c < nMultiConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH1F(Form("h1d_%s_layer%d", multi_categories[c], layer),
                            Form("Layer %d: %s;ToT [ns];Hits", layer, multi_categories[c]),
                            nBins, xMin, xMax);
            strip_histograms[multi_categories[c]][layer] = hist;
        }
    }

    // 2D Correlation Histograms (Primary ToT vs Secondary ToT)
    std::map<std::string, std::map<int, TH2F*>> h2d_tot_corr;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            h2d_tot_corr[categories[c]][layer] = new TH2F(
                Form(
                    "h2d_corr_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Primary ToT [ns];Secondary ToT [ns];Entries", layer, comments[c]),
                nBins, xMin, xMax, nBins, xMin, xMax);
        }
    }

    while (readerProcData.Next() && readerTrackData.Next()) {
        const size_t n_hits = std::min({tot1->size(), tot2->size(), strips->size(), layers->size(), dts->size(), time1->size(), time2->size()});

        // First pass: Find minimum processed hit time per strip per side for the current event
        std::map<std::pair<int, int>, int> min_time_eta1;
        std::map<std::pair<int, int>, int> min_time_eta2;

        for (size_t i = 0; i < n_hits; ++i) {
            int layer = (*layers)[i];
            if (layer < 0 || layer > 2) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            auto key = std::make_pair(layer, strip);

            if ((*has_valid_time_eta1)[i]) {
                int t1 = (*time1)[i];
                if (min_time_eta1.find(key) == min_time_eta1.end() || t1 < min_time_eta1[key]) {
                    min_time_eta1[key] = t1;
                }
            }
            if ((*has_valid_time_eta2)[i]) {
                int t2 = (*time2)[i];
                if (min_time_eta2.find(key) == min_time_eta2.end() || t2 < min_time_eta2[key]) {
                    min_time_eta2[key] = t2;
                }
            }
        }

        // Pass 1.5: Extract Primary ToT for 2D correlations
        std::map<std::pair<int, int>, double> primary_tot_eta1;
        std::map<std::pair<int, int>, double> primary_tot_eta2;

        for (size_t i = 0; i < n_hits; ++i) {
            int layer = (*layers)[i];
            if (layer < 0 || layer > 2) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            auto key = std::make_pair(layer, strip);

            if ((*has_valid_time_eta1)[i] && (*time1)[i] == min_time_eta1[key]) {
                primary_tot_eta1[key] = TimeUtils::ticksToTime((*tot1)[i]);
            }
            if ((*has_valid_time_eta2)[i] && (*time2)[i] == min_time_eta2[key]) {
                primary_tot_eta2[key] = TimeUtils::ticksToTime((*tot2)[i]);
            }
        }

        // Second pass: Fill 1D histograms and 2D correlation histograms
        for (size_t i = 0; i < n_hits; ++i) {
            int layer = (*layers)[i];
            if (layer < 0 || layer > 2) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            auto key = std::make_pair(layer, strip);
            int dt = (*dts)[i];

            double tot1_ns = TimeUtils::ticksToTime((*tot1)[i]);
            double tot2_ns = TimeUtils::ticksToTime((*tot2)[i]);
            bool valid_track = (*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i];

            bool in_beam = std::find(region_of_interest.beam_region.begin(),
                                     region_of_interest.beam_region.end(),
                                     std::make_pair(strip, dt)) != region_of_interest.beam_region.end();

            // Side Eta1
            if ((*has_valid_time_eta1)[i]) {
                bool is_primary_eta1 = ((*time1)[i] == min_time_eta1[key]);
                std::string tag_eta1 = is_primary_eta1 ? "primary" : "secondary";

                strip_histograms["tot_eta1_before_reco"][layer]->Fill(tot1_ns);
                strip_histograms["tot_eta1_" + tag_eta1 + "_before_reco"][layer]->Fill(tot1_ns);

                if (!is_primary_eta1 && primary_tot_eta1.count(key)) {
                    double prim_tot = primary_tot_eta1[key];
                    h2d_tot_corr["tot_eta1_before_reco"][layer]->Fill(prim_tot, tot1_ns);
                    if (valid_track) h2d_tot_corr["tot_eta1_after_reco"][layer]->Fill(prim_tot, tot1_ns);
                    if (in_beam) h2d_tot_corr["tot_eta1_beam"][layer]->Fill(prim_tot, tot1_ns);
                }

                if (valid_track) {
                    strip_histograms["tot_eta1_after_reco"][layer]->Fill(tot1_ns);
                    strip_histograms["tot_eta1_" + tag_eta1 + "_after_reco"][layer]->Fill(tot1_ns);
                }
                if (in_beam) {
                    strip_histograms["tot_eta1_beam"][layer]->Fill(tot1_ns);
                    strip_histograms["tot_eta1_" + tag_eta1 + "_beam"][layer]->Fill(tot1_ns);
                }
            }

            // Side Eta2
            if ((*has_valid_time_eta2)[i]) {
                bool is_primary_eta2 = ((*time2)[i] == min_time_eta2[key]);
                std::string tag_eta2 = is_primary_eta2 ? "primary" : "secondary";

                strip_histograms["tot_eta2_before_reco"][layer]->Fill(tot2_ns);
                strip_histograms["tot_eta2_" + tag_eta2 + "_before_reco"][layer]->Fill(tot2_ns);

                if (!is_primary_eta2 && primary_tot_eta2.count(key)) {
                    double prim_tot = primary_tot_eta2[key];
                    h2d_tot_corr["tot_eta2_before_reco"][layer]->Fill(prim_tot, tot2_ns);
                    if (valid_track) h2d_tot_corr["tot_eta2_after_reco"][layer]->Fill(prim_tot, tot2_ns);
                    if (in_beam) h2d_tot_corr["tot_eta2_beam"][layer]->Fill(prim_tot, tot2_ns);
                }

                if (valid_track) {
                    strip_histograms["tot_eta2_after_reco"][layer]->Fill(tot2_ns);
                    strip_histograms["tot_eta2_" + tag_eta2 + "_after_reco"][layer]->Fill(tot2_ns);
                }
                if (in_beam) {
                    strip_histograms["tot_eta2_beam"][layer]->Fill(tot2_ns);
                    strip_histograms["tot_eta2_" + tag_eta2 + "_beam"][layer]->Fill(tot2_ns);
                }
            }
        }
    }

    double global_percentile = 0.0;
    double percentile = 0.99;

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            TH1* h = strip_histograms[categories[c]][layer];
            if (h->GetEffectiveEntries() > 0) {
                double q[1];
                double prob[1] = {percentile};
                h->GetQuantiles(1, q, prob);
                if (q[0] > global_percentile) {
                    global_percentile = q[0];
                }
            }
        }
    }

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            strip_histograms[categories[c]][layer]->GetXaxis()->SetRangeUser(0.0, global_percentile);
            h2d_tot_corr[categories[c]][layer]->GetXaxis()->SetRangeUser(0.0, global_percentile);
            h2d_tot_corr[categories[c]][layer]->GetYaxis()->SetRangeUser(0.0, global_percentile);
        }
    }
    for (int c = 0; c < nMultiConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            strip_histograms[multi_categories[c]][layer]->GetXaxis()->SetRangeUser(0.0, global_percentile);
        }
    }

    struct StackPairing {
        std::string eta1_cat;
        std::string eta2_cat;
        std::string suffix;
        std::string title_modifier;
    };

    std::vector<StackPairing> pairings = {
        {"tot_eta1_before_reco", "tot_eta2_before_reco", "before_reco", "Before Track Reco"},
        {"tot_eta1_after_reco",  "tot_eta2_after_reco",  "after_reco",  "After Track Reco"},
        {"tot_eta1_beam",        "tot_eta2_beam",        "beam",        "Beam Spot Region"}
    };

    // Store stacks combining both sides (eta1 vs eta2 per layer)
    for (int layer : {0, 1, 2}) {
        for (const auto& pair : pairings) {
            auto* h_tot1 = strip_histograms[pair.eta1_cat][layer];
            auto* h_tot2 = strip_histograms[pair.eta2_cat][layer];

            h_tot1->SetLineColor(kBlue);
            h_tot1->SetMarkerColor(kBlue);
            
            h_tot2->SetLineColor(kRed);
            h_tot2->SetMarkerColor(kRed);

            std::string stack_name = Form(
                "stack_tot_side_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                pair.suffix.c_str(),
                layer,
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            );
            std::string stack_title = Form("Layer %d: %s", layer, pair.title_modifier.c_str());

            auto* stack = new THStack(stack_name.c_str(), stack_title.c_str());

            stack->Add(h_tot2);
            stack->Add(h_tot1);

            stack->Write("", TObject::kOverwrite);
            delete stack;
        }
    }

    // Multi-hit comparison stacks per side (primary + secondary hits)
    for (int layer : {0, 1, 2}) {
        for (const auto& pair : pairings) {
            // Side Eta1 Multi-hit Stack
            auto* h_e1_primary   = strip_histograms["tot_eta1_primary_" + pair.suffix][layer];
            auto* h_e1_secondary = strip_histograms["tot_eta1_secondary_" + pair.suffix][layer];

            std::string stack_e1_name = Form(
                "stack_tot_multi_eta1_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                pair.suffix.c_str(),
                layer,
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            );
            std::string stack_e1_title = Form("Layer %d: Side #eta+: %s", layer, pair.title_modifier.c_str());

            auto* stack_e1 = new THStack(stack_e1_name.c_str(), stack_e1_title.c_str());
            stack_e1->Add(h_e1_primary);
            stack_e1->Add(h_e1_secondary);
            stack_e1->Write("", TObject::kOverwrite);
            delete stack_e1;

            // Side Eta2 Multi-hit Stack
            auto* h_e2_primary   = strip_histograms["tot_eta2_primary_" + pair.suffix][layer];
            auto* h_e2_secondary = strip_histograms["tot_eta2_secondary_" + pair.suffix][layer];

            std::string stack_e2_name = Form(
                "stack_tot_multi_eta2_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                pair.suffix.c_str(),
                layer,
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            );
            std::string stack_e2_title = Form("Layer %d: Side #eta-: %s", layer, pair.title_modifier.c_str());

            auto* stack_e2 = new THStack(stack_e2_name.c_str(), stack_e2_title.c_str());
            stack_e2->Add(h_e2_primary);
            stack_e2->Add(h_e2_secondary);
            stack_e2->Write("", TObject::kOverwrite);
            delete stack_e2;
        }
    }

    // Store stacks combining layers
    for (int c = 0; c < nConfigs; ++c) {
        std::string cat = categories[c];
        std::string side_label = (cat.find("eta1") != std::string::npos) ? "Side #eta+" : "Side #eta-";
        std::string stack_name = Form(
            "stack_tot_layers_%s_%s_lv_%d_filter_%s_hv_%.0f",
            cat.c_str(),
            _metadata.mixture.c_str(),
            _metadata.lv_setting,
            filter_str.c_str(),
            _metadata.scanned_hv
        );
        std::string stack_title = Form("%s: %s", comments[c], side_label.c_str());

        auto* stack = new THStack(stack_name.c_str(), stack_title.c_str());

        for (int layer : {0, 1, 2}) {
            stack->Add(strip_histograms[cat][layer]);
        }

        stack->Write("", TObject::kOverwrite);
        delete stack;
    }

    // Write all individual 1D histograms
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            strip_histograms[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete strip_histograms[categories[c]][layer];
        }
    }
    for (int c = 0; c < nMultiConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            strip_histograms[multi_categories[c]][layer]->Write("", TObject::kOverwrite);
            delete strip_histograms[multi_categories[c]][layer];
        }
    }

    // Write and clean up 2D correlation histograms
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            h2d_tot_corr[categories[c]][layer]->Write("", TObject::kOverwrite);
            delete h2d_tot_corr[categories[c]][layer];
        }
    }
}

void plotToTVsStrip(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* tot_strip_dir = analysis_dir->GetDirectory("tot_strip");
    tot_strip_dir->cd();

    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!proc_tree || proc_tree->IsZombie() || !track_tree || track_tree->IsZombie()) {
        std::cerr << "Error: Invalid trees for plotToTVsStrip analysis." << std::endl;
        return;
    }

    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> proc_tot1(readerProcData, "proc_tot1");
    TTreeReaderValue<std::vector<int>> proc_tot2(readerProcData, "proc_tot2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    const int nConfigs = 4;
    const char* suffixes[nConfigs] = {
        "tot_eta1_before_reco", "tot_eta2_before_reco", 
        "tot_eta1_after_reco", "tot_eta2_after_reco"
    };
    const char* comments[nConfigs] = {
        "Side #eta+: Before Track Reco", "Side #eta-: Before Track Reco",
        "Side #eta+: After Track Reco", "Side #eta-: After Track Reco"
    };

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    std::map<std::string, std::map<int, TH2*>> tot_strip_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TH2F(
                Form("h2d_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    suffixes[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Strip;ToT [ns];Entries", layer, comments[c]),
                STRIPS_PER_LAYER, 0, STRIPS_PER_LAYER, 25, 0, 25
            );
            tot_strip_histograms[suffixes[c]][layer] = hist;
        }
    }

    while (readerProcData.Next() && readerTrackData.Next()) {
        const size_t n_hits = std::min({strips->size(), layers->size(), proc_tot1->size(), proc_tot2->size()});
        for (size_t i = 0; i < n_hits; ++i) {
            int layer = (*layers)[i];
            if (layer < 0 || layer > 2) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            double tot1_ns = TimeUtils::ticksToTime((*proc_tot1)[i]);
            double tot2_ns = TimeUtils::ticksToTime((*proc_tot2)[i]);
            bool valid_track = (*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i];

            // Before Reco
            tot_strip_histograms["tot_eta1_before_reco"][layer]->Fill(strip, tot1_ns);
            tot_strip_histograms["tot_eta2_before_reco"][layer]->Fill(strip, tot2_ns);

            // After Reco (either eta1 or eta2 valid track)
            if (valid_track) {
                tot_strip_histograms["tot_eta1_after_reco"][layer]->Fill(strip, tot1_ns);
                tot_strip_histograms["tot_eta2_after_reco"][layer]->Fill(strip, tot2_ns);
            }
        }
    }

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {   
            tot_strip_histograms[suffixes[c]][layer]->Write("", TObject::kOverwrite);
            delete tot_strip_histograms[suffixes[c]][layer];
        }
    }
}

void plotToTVsDtVsStrip(TFile* input_file, const RoI& region_of_interest) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* tot_dt_strip_dir = analysis_dir->GetDirectory("tot_dt_strip");
    if (!tot_dt_strip_dir) {
        tot_dt_strip_dir = analysis_dir->mkdir("tot_dt_strip");
    }
    tot_dt_strip_dir->cd();

    TTree* input_data_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!input_data_tree || !proc_tree || !track_tree) {
        std::cerr << "Error: Invalid trees for plotToTVsDtVsStrip analysis." << std::endl;
        return;
    }

    TTreeReader readerInputData(input_data_tree);
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);

    TTreeReaderValue<std::vector<int>> raw_time1(readerInputData, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>> raw_time2(readerInputData, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> dts(readerProcData, "proc_dt_time1_time2");
    TTreeReaderValue<std::vector<int>> proc_tot1(readerProcData, "proc_tot1");
    TTreeReaderValue<std::vector<int>> proc_tot2(readerProcData, "proc_tot2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    const int nConfigs = 6;
    const char* categories[nConfigs] = {
        "tot_eta1_dt_strip_before_reco", "tot_eta2_dt_strip_before_reco",
        "tot_eta1_dt_strip_after_reco", "tot_eta2_dt_strip_after_reco",
        "tot_eta1_dt_strip_beam", "tot_eta2_dt_strip_beam"
    };
    const char* comments[nConfigs] = {
        "Side #eta+: Before Track Reco", "Side #eta-: Before Track Reco",
        "Side #eta+: After Track Reco", "Side #eta-: After Track Reco",
        "Side #eta+: Beam Spot Region", "Side #eta-: Beam Spot Region"
    };

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    std::map<std::string, std::map<int, TProfile2D*>> profile_histograms;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            auto* hist = new TProfile2D(
                Form(
                    "h2d_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv),
                Form("Layer %d: %s;Strip;#Delta#it{t} [ns];<ToT> [ns]", layer, comments[c]),
                STRIPS_PER_LAYER, 0, STRIPS_PER_LAYER, 20, -10, 10
            );
            profile_histograms[categories[c]][layer] = hist;
        }
    }

    while (readerInputData.Next() && readerProcData.Next() && readerTrackData.Next()) {
        const size_t n_hits = std::min({strips->size(), layers->size(), dts->size(), proc_tot1->size(), proc_tot2->size()});
        for (size_t i = 0; i < n_hits; ++i) {
            int layer = (*layers)[i];
            if (layer < 0 || layer > 2) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            int dt = TimeUtils::ticksToTime((*dts)[i]);
            double tot1_ns = TimeUtils::ticksToTime((*proc_tot1)[i]);
            double tot2_ns = TimeUtils::ticksToTime((*proc_tot2)[i]);
            bool valid_track = (*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i];

            bool in_beam = std::find(region_of_interest.beam_region.begin(),
                                     region_of_interest.beam_region.end(),
                                     std::make_pair(strip, dt)) != region_of_interest.beam_region.end();

            // Before Reco
            if ((*raw_time1)[i] != 0) profile_histograms["tot_eta1_dt_strip_before_reco"][layer]->Fill(strip, dt, tot1_ns);
            if ((*raw_time2)[i] != 0) profile_histograms["tot_eta2_dt_strip_before_reco"][layer]->Fill(strip, dt, tot2_ns);

            // After Reco
            if ((*raw_time1)[i] != 0 && valid_track) profile_histograms["tot_eta1_dt_strip_after_reco"][layer]->Fill(strip, dt, tot1_ns);
            if ((*raw_time2)[i] != 0 && valid_track) profile_histograms["tot_eta2_dt_strip_after_reco"][layer]->Fill(strip, dt, tot2_ns);

            // Beam Region
            if ((*raw_time1)[i] != 0 && in_beam) profile_histograms["tot_eta1_dt_strip_beam"][layer]->Fill(strip, dt, tot1_ns);
            if ((*raw_time2)[i] != 0 && in_beam) profile_histograms["tot_eta2_dt_strip_beam"][layer]->Fill(strip, dt, tot2_ns);
        }
    }

    const double min_entries_cut = 5.0;
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            TProfile2D* prof = profile_histograms[categories[c]][layer];

            int n_bins_x = prof->GetNbinsX();
            int n_bins_y = prof->GetNbinsY();

            for (int x = 1; x <= n_bins_x; ++x) {
                for (int y = 1; y <= n_bins_y; ++y) {
                    int bin = prof->GetBin(x, y);
                    if (prof->GetBinEntries(bin) < min_entries_cut) {
                        prof->SetBinContent(bin, 0.0);
                        prof->SetBinEntries(bin, 0.0);
                    }
                }
            }

            prof->SetOption("COLZ");
            prof->Write("", TObject::kOverwrite);
            delete prof;
        }
    }
}

void plotMultiplicity(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* mult_strip_dir = analysis_dir->GetDirectory("multiplicity_strip");
    if (!mult_strip_dir) {
        mult_strip_dir = analysis_dir->mkdir("multiplicity_strip");
    }
    mult_strip_dir->cd();

    TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));

    if (!input_tree || !proc_tree || !track_tree) {
        std::cerr << "Error: Invalid trees for plotMultiplicity analysis." << std::endl;
        return;
    }

    TTreeReader readerInputData(input_tree);
    TTreeReader readerProcData(proc_tree);
    TTreeReader readerTrackData(track_tree);
    TTreeReaderValue<std::vector<int>> raw_time1(readerInputData, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>> raw_time2(readerInputData, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>> rise(readerInputData, "hit_rise");
    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    const int nConfigs = 4;
    const char* categories[nConfigs] = {
        "mult_eta1_before_reco", "mult_eta2_before_reco",
        "mult_eta1_after_reco",  "mult_eta2_after_reco"
    };
    const char* comments[nConfigs] = {
        "Side #eta+: Before Track Reco", "Side #eta-: Before Track Reco",
        "Side #eta+: After Track Reco", "Side #eta-: After Track Reco"
    };

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    // [category][layer] -> histogram
    std::map<std::string, std::map<int, TH2F*>> h2d_mult;         // Strip vs Multiplicity vs Entries
    std::map<std::string, std::map<int, TH1F*>> h1d_mult;         // Multiplicity vs Entries (Integrated)
    std::map<std::string, std::map<int, TProfile*>> prof_avg;     // Strip vs Average Multiplicity
    std::map<std::string, std::map<int, TProfile*>> prof_frac;    // Strip vs Fraction of Mult > 1

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            // 2D Map (Strip vs Multiplicity)
            h2d_mult[categories[c]][layer] = new TH2F(
                Form(
                    "h2d_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Strip;Multiplicity;Entries", layer, comments[c]),
                STRIPS_PER_LAYER, 0, STRIPS_PER_LAYER, 9, 1, 10);

            // 1D Distribution (Summed Multiplicity)
            h1d_mult[categories[c]][layer] = new TH1F(
                Form(
                    "h1d_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Multiplicity;Entries", layer, comments[c]),
                4, 1, 5);

            // 1D Profile (Average Multiplicity vs Strip)
            prof_avg[categories[c]][layer] = new TProfile(
                Form(
                    "prof_avg_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Strip;#LTAverage Multiplicity#GT [Hits]", layer, comments[c]),
                STRIPS_PER_LAYER, 0, STRIPS_PER_LAYER);

            // 1D Profile (Fraction of Mult > 1 vs Strip)
            prof_frac[categories[c]][layer] = new TProfile(
                Form(
                    "prof_frac_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Strip;#frac{N_{mult > 1}}{N_{events}}", layer, comments[c]),
                STRIPS_PER_LAYER, 0, STRIPS_PER_LAYER);
        }
    }

    // Event Loop
    while (readerInputData.Next() && readerProcData.Next() && readerTrackData.Next()) {

        // Reset lookup table for the current event: counts[category][layer][strip]
        std::map<int, std::map<int, int>> counts[nConfigs];

        const size_t n_hits = std::min({
            strips->size(),
            layers->size(),
            rise->size(),
            raw_time1->size(),
            raw_time2->size(),
            in_valid_track_eta1->size(),
            in_valid_track_eta2->size()
        });

        // Count hits per strip for this event
        for (size_t i = 0; i < n_hits; ++i) {
            if ((*rise)[i] == 0) continue; // Skip falling edge hits

            int layer = (*layers)[i];
            if (layer < 0 || layer > 2) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            bool valid_track = (*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i];

            if ((*raw_time1)[i] != 0) {
                counts[0][layer][strip]++;
                if (valid_track) counts[2][layer][strip]++;
            }
            if ((*raw_time2)[i] != 0) {
                counts[1][layer][strip]++;
                if (valid_track) counts[3][layer][strip]++;
            }
        }

        // Fill all 4 statistics using the aggregated counts
        for (int c = 0; c < nConfigs; ++c) {
            for (int layer : {0, 1, 2}) {
                for (const auto& [strip, count] : counts[c][layer]) {
                    // Only process strips that actually had hits this event
                    if (count > 0) {
                        h2d_mult[categories[c]][layer]->Fill(strip, count);
                        h1d_mult[categories[c]][layer]->Fill(count);
                        prof_avg[categories[c]][layer]->Fill(strip, count);

                        // Pass 1.0 if mult > 1, pass 0.0 if mult == 1.
                        prof_frac[categories[c]][layer]->Fill(strip, count > 1 ? 1.0 : 0.0);
                    }
                }
            }
        }
    }

    mult_strip_dir->cd();

    int colors[3] = {kBlack, kBlue+1, kRed+1};
    int markers[3] = {20, 21, 22};

    for (int c = 0; c < nConfigs; ++c) {
        // Write the 2D maps and 1D distributions individually
        for (int layer : {0, 1, 2}) {
            h2d_mult[categories[c]][layer]->Write("", TObject::kOverwrite);
            h1d_mult[categories[c]][layer]->Write("", TObject::kOverwrite);
        }

        // Combine the Average Multiplicity Profiles into a THStack
        THStack* stack_avg = new THStack(
            Form(
                "stack_avg_%s_%s_lv_%d_filter_%s_hv_%.0f",
                categories[c],
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            ),
            Form("%s;Strip;Average Multiplicity", comments[c])
        );
        for (int layer : {0, 1, 2}) {
            prof_avg[categories[c]][layer]->SetLineColor(colors[layer]);
            prof_avg[categories[c]][layer]->SetMarkerColor(colors[layer]);
            prof_avg[categories[c]][layer]->SetMarkerStyle(markers[layer]);

            stack_avg->Add(prof_avg[categories[c]][layer], "PE");
        }
        stack_avg->Write("", TObject::kOverwrite);

        // Combine the Fraction Profiles into a THStack
        THStack* stack_frac = new THStack(
            Form(
                "stack_frac_%s_%s_lv_%d_filter_%s_hv_%.0f",
                categories[c],
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            ),
            Form("%s;Strip;N_{mult > 1} / N_{events}", comments[c])
        );
        for (int layer : {0, 1, 2}) {
            prof_frac[categories[c]][layer]->SetLineColor(colors[layer]);
            prof_frac[categories[c]][layer]->SetMarkerColor(colors[layer]);
            prof_frac[categories[c]][layer]->SetMarkerStyle(markers[layer]);

            stack_frac->Add(prof_frac[categories[c]][layer], "PE");
        }
        stack_frac->Write("", TObject::kOverwrite);

        delete stack_avg;
        delete stack_frac;

        for (int layer : {0, 1, 2}) {
            delete h2d_mult[categories[c]][layer];
            delete h1d_mult[categories[c]][layer];
            delete prof_avg[categories[c]][layer];
            delete prof_frac[categories[c]][layer];
        }
    }
}

void plotDelay(TFile* input_file) {
    TDirectory* analysis_dir = input_file->GetDirectory("analysis");
    TDirectory* delay_strip_dir = analysis_dir->GetDirectory("delay_strip");
    if (!delay_strip_dir) {
        delay_strip_dir = analysis_dir->mkdir("delay_strip");
    }
    delay_strip_dir->cd();

    TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("InputData"));
    TTree* proc_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));

    if (!input_tree || !proc_tree || !track_tree) {
        std::cerr << "Error: Invalid trees for plotDelay analysis." << std::endl;
        return;
    }

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

    const int nConfigs = 4;
    const char* categories[nConfigs] = {
        "delay_eta1_before_reco", "delay_eta2_before_reco",
        "delay_eta1_after_reco",  "delay_eta2_after_reco"
    };
    const char* comments[nConfigs] = {
        "Side #eta+: Before Track Reco", "Side #eta-: Before Track Reco",
        "Side #eta+: After Track Reco", "Side #eta-: After Track Reco"
    };

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    std::map<std::string, std::map<int, TH2F*>> h2d_delay;
    std::map<std::string, std::map<int, std::vector<TH1F*>>> h1d_delay;
    std::map<std::string, std::map<int, TProfile*>> prof_delay;

    const int max_delay_bins = 100;
    double x_min = -0.5 * TIME_TICK_NS;
    double x_max = (max_delay_bins - 0.5) * TIME_TICK_NS;

    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            h2d_delay[categories[c]][layer] = new TH2F(
                Form(
                    "h2d_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Strip;Delay from First Hit [ns];Entries", layer, comments[c]),
                STRIPS_PER_LAYER, 0, STRIPS_PER_LAYER, max_delay_bins, x_min, x_max);

            prof_delay[categories[c]][layer] = new TProfile(
                Form(
                    "prof_avg_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Strip;#LTAverage Delay#GT [ns]", layer, comments[c]),
                STRIPS_PER_LAYER, 0, STRIPS_PER_LAYER);

            // Create a 3-element vector of 1D Histograms for 2nd, 3rd, and 4th+ hits
            for (int hit = 0; hit < 3; ++hit) {
                std::string hit_label = (hit == 0) ? "2nd Hit" : (hit == 1) ? "3rd Hit" : "4th+ Hit";

                h1d_delay[categories[c]][layer].push_back(new TH1F(
                    Form("h1d_%s_layer%d_hit%d", categories[c], layer, hit + 2),
                    hit_label.c_str(),
                    max_delay_bins, x_min, x_max));
            }
        }
    }

    while (readerInputData.Next() && readerProcData.Next() && readerTrackData.Next()) {
        std::map<int, std::map<int, std::vector<int>>> event_delays[nConfigs];

        const size_t n_hits = std::min({strips->size(), layers->size(), rise->size(), time1->size(), time2->size()});

        for (size_t i = 0; i < n_hits; ++i) {
            if ((*rise)[i] == 0) continue;

            int layer = (*layers)[i];
            if (layer < 0 || layer > 2) continue;

            int strip = perFileHelpers::remapStrip((*strips)[i]);
            bool valid_track = (*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i];

            if ((*raw_time1)[i] != 0) {
                event_delays[0][layer][strip].push_back((*time1)[i]);
                if (valid_track) event_delays[2][layer][strip].push_back((*time1)[i]);
            }
            if ((*raw_time2)[i] != 0) {
                event_delays[1][layer][strip].push_back((*time2)[i]);
                if (valid_track) event_delays[3][layer][strip].push_back((*time2)[i]);
            }
        }

        for (int c = 0; c < nConfigs; ++c) {
            for (int layer : {0, 1, 2}) {
                for (const auto& [strip, delay_vec_raw] : event_delays[c][layer]) {

                    if (delay_vec_raw.size() > 1) {
                        std::vector<int> sorted_delays = delay_vec_raw;
                        std::sort(sorted_delays.begin(), sorted_delays.end());

                        int min_delay = sorted_delays[0];

                        for (size_t i = 1; i < sorted_delays.size(); ++i) {
                            double dt = TimeUtils::ticksToTime(sorted_delays[i] - min_delay);

                            h2d_delay[categories[c]][layer]->Fill(strip, dt);
                            prof_delay[categories[c]][layer]->Fill(strip, dt);

                            int hit_idx = std::min(static_cast<int>(i) - 1, 2);
                            h1d_delay[categories[c]][layer][hit_idx]->Fill(dt);
                        }
                    }
                }
            }
        }
    }

    double percentile = 0.99;
    double global_percentile = 0.0;

    // Find the 99th percentile across all populated 1D histograms
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            for (int hit = 0; hit < 3; ++hit) {
                TH1F* h = h1d_delay[categories[c]][layer][hit];
                if (h->GetEffectiveEntries() > 0) {
                    double q[1];
                    double prob[1] = {percentile};
                    h->GetQuantiles(1, q, prob);
                    if (q[0] > global_percentile) {
                        global_percentile = q[0];
                    }
                }
            }
        }
    }

    // Add visual buffer (2 ticks worth of ns) and clamp to absolute axis max
    double dynamic_max = std::min(global_percentile + (2.0 * TIME_TICK_NS), x_max);

    // Apply the limit to the raw plots before writing
    for (int c = 0; c < nConfigs; ++c) {
        for (int layer : {0, 1, 2}) {
            // Apply to 2D Heatmap (Y-Axis)
            h2d_delay[categories[c]][layer]->GetYaxis()->SetRangeUser(0.0, dynamic_max);

            // Apply to 1D Distributions (X-Axis)
            for (int hit = 0; hit < 3; ++hit) {
                h1d_delay[categories[c]][layer][hit]->GetXaxis()->SetRangeUser(0.0, dynamic_max);
            }
        }
    }

    delay_strip_dir->cd();

    int colors[3] = {kBlack, kBlue+1, kRed+1};
    int markers[3] = {20, 21, 22};

    for (int c = 0; c < nConfigs; ++c) {
        THStack* stack_avg = new THStack(
            Form(
                "stack_avg_%s_%s_lv_%d_filter_%s_hv_%.0f",
                categories[c],
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            ),
                Form("%s;Strip;Average Delay [ns]", comments[c])
            );

        for (int layer : {0, 1, 2}) {
            h2d_delay[categories[c]][layer]->Write("", TObject::kOverwrite);

            // Package the 1D Distributions into a THStack
            THStack* stack_dist = new THStack(
                Form(
                    "stack_%s_layer%d_%s_lv_%d_filter_%s_hv_%.0f",
                    categories[c],
                    layer,
                    _metadata.mixture.c_str(),
                    _metadata.lv_setting,
                    filter_str.c_str(),
                    _metadata.scanned_hv
                ),
                Form("Layer %d: %s;Delay from First Hit [ns];Entries", layer, comments[c])
            );

            for (int hit = 0; hit < 3; ++hit) {
                stack_dist->Add(h1d_delay[categories[c]][layer][hit]);
            }
            stack_dist->Write("", TObject::kOverwrite);
            delete stack_dist;

            prof_delay[categories[c]][layer]->SetLineColor(colors[layer]);
            prof_delay[categories[c]][layer]->SetMarkerColor(colors[layer]);
            prof_delay[categories[c]][layer]->SetMarkerStyle(markers[layer]);

            stack_avg->Add(prof_delay[categories[c]][layer], "PE");
        }

        stack_avg->Write("", TObject::kOverwrite);

        delete stack_avg;
        for (int layer : {0, 1, 2}) {
            delete h2d_delay[categories[c]][layer];
            for (int hit = 0; hit < 3; ++hit) {
                delete h1d_delay[categories[c]][layer][hit];
            }
            delete prof_delay[categories[c]][layer];
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

    // Buffers to accumulate raw ToF values (in ns) prior to histogram creation
    std::array<std::vector<double>, LAYER_PAIR_COUNT> raw_tof1_data;
    std::array<std::vector<double>, LAYER_PAIR_COUNT> raw_tof2_data;

    std::string filter_str = Form("%.1f", _metadata.filter);
    std::replace(filter_str.begin(), filter_str.end(), '.', '_');

    for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
        std::string branch1 = "track_time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta1";
        std::string branch2 = "track_time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i] + "_eta2";

        tof1_readers[i] = std::make_unique<TTreeReaderValue<std::vector<int>>>(readerTrackData, branch1.c_str());
        tof2_readers[i] = std::make_unique<TTreeReaderValue<std::vector<int>>>(readerTrackData, branch2.c_str());
    }

    while (readerTrackData.Next()) {
        for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
            if (tof1_readers[i]->GetSetupStatus() == 0) {
                for (int t : **tof1_readers[i]) raw_tof1_data[i].push_back(TimeUtils::ticksToTime(t));
            }
            if (tof2_readers[i]->GetSetupStatus() == 0) {
                for (int t : **tof2_readers[i]) raw_tof2_data[i].push_back(TimeUtils::ticksToTime(t));
            }
        }
    }

    // Helper lambda to apply the visual fit before saving
    auto applyVisualFit = [](TH1F* hist) {
        if (!hist || hist->Integral() < 20) return;

        TFitResultPtr r1 = hist->Fit("gaus", "QS");

        if (r1.Get() != nullptr && r1->IsValid() && static_cast<int>(r1) == 0) {

            double chi2 = r1->Chi2();
            double ndf  = r1->Ndf();
            double chi2_ndf = (ndf > 0) ? (chi2 / ndf) : 999.0;

            if (chi2_ndf < 2.5) return;

            double p_mean  = r1->Parameter(1);
            double p_sigma = r1->Parameter(2);

            if (p_sigma > 0.0) {
                double fit_width = std::max(1.5 * p_sigma, 2.0 * TIME_TICK_NS);

                TFitResultPtr r2 = hist->Fit("gaus", "QS", "", p_mean - fit_width, p_mean + fit_width);

                if (r2.Get() == nullptr || !r2->IsValid() || static_cast<int>(r2) != 0) {
                    hist->Fit("gaus", "QS");
                }

                if (TF1* fit_func = hist->GetFunction("gaus")) {
                    fit_func->SetRange(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
                }
            }
        }
    };

    double bin_width = TIME_TICK_NS;
    for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
        auto [x_min1, x_max1] = perFileHelpers::getPercentileRange(raw_tof1_data[i], 0.0, 1.0);
        auto [x_min2, x_max2] = perFileHelpers::getPercentileRange(raw_tof2_data[i], 0.0, 1.0);

        int n_bins1 = std::max(1, static_cast<int>(std::round((x_max1 - x_min1) / bin_width)));
        int n_bins2 = std::max(1, static_cast<int>(std::round((x_max2 - x_min2) / bin_width)));

        x_max1 = x_min1 + n_bins1 * bin_width;
        x_max2 = x_min2 + n_bins2 * bin_width;

        TH1F* h_tof1 = new TH1F(
            Form(
                "h1d_tof_layer_pair_%s_layer%d_eta1_%s_lv_%d_filter_%s_hv_%.0f",
                LAYER_PAIR_SUFFIXES[i].c_str(),
                _metadata.scanned_layer,
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            ),
            Form("Side #eta+: %s;ToF [ns];Entries", pair_labels[i].c_str()),
            n_bins1, x_min1, x_max1
        );

        TH1F* h_tof2 = new TH1F(
            Form(
                "h1d_tof_layer_pair_%s_layer%d_eta2_%s_lv_%d_filter_%s_hv_%.0f",
                LAYER_PAIR_SUFFIXES[i].c_str(),
                _metadata.scanned_layer,
                _metadata.mixture.c_str(),
                _metadata.lv_setting,
                filter_str.c_str(),
                _metadata.scanned_hv
            ),
            Form("Side #eta-: %s;ToF [ns];Entries", pair_labels[i].c_str()),
            n_bins2, x_min2, x_max2
        );

        for (double val : raw_tof1_data[i]) h_tof1->Fill(val);
        for (double val : raw_tof2_data[i]) h_tof2->Fill(val);

        applyVisualFit(h_tof1);
        applyVisualFit(h_tof2);

        h_tof1->Write("", TObject::kOverwrite);
        h_tof2->Write("", TObject::kOverwrite);

        delete h_tof1;
        delete h_tof2;
    }
}

}