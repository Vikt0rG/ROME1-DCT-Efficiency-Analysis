#include "process_raw_hits.hpp"
#include <TROOT.h>
#include <iostream>
#include <cmath>

// ============================================================================================================
// Setup functions

void setupOutputFile(TFile*& output_file, 
                     TTree*& input_data_tree, 
                     TTree*& processed_data_tree, 
                     TTree*& clusterization_tree) {
    // Create output ROOT file
    output_file = TFile::Open("out.root", "RECREATE");
    if (!output_file) {
        std::cerr << "ERROR: Failed to create output ROOT file!" << std::endl;
        return;
    }
    
    // Create trees
    input_data_tree = new TTree("input_data", "Input DCT data");
    processed_data_tree = new TTree("processed_data", "Processed DCT data");
    clusterization_tree = new TTree("clusterization_data", "Data after clusterization and corrections");
}

void initializeEfficiencyCounters(EfficiencyCounters& counters) {
    counters.triggered_events_external = 0;
    
    for (int i = 0; i < 3; i++) {
        counters.triggered_events_rpc[i] = 0;
        counters.eta1_efficiency_counter[i] = 0;
        counters.eta2_efficiency_counter[i] = 0;
        counters.eta_or_efficiency_counter[i] = 0;
        counters.eta_and_efficiency_counter[i] = 0;
        counters.eta1_efficiency_counter_rpc[i] = 0;
        counters.eta2_efficiency_counter_rpc[i] = 0;
        counters.eta_or_efficiency_counter_rpc[i] = 0;
        counters.eta_and_efficiency_counter_rpc[i] = 0;
    }
}

// ============================================================================================================
// Word processing functions

void extractWordData(int word, int& channel, int& rise, int& raw_bcid, int& raw_time1, int& raw_time2) {
    channel = (word >> 20) & 0xFF;    // Bits 20-27
    rise = word & 0x01;               // Bit 0

    if (rise == 1) {                  // Rising edge
        raw_bcid = word >> 12 & 0xFF;   // Bits 12-19
        raw_time1 = word >> 7 & 0x1F;   // Bits 7-11
        raw_time2 = word >> 1 & 0x3F;   // Bits 1-5
    } else {                          // Falling edge
        raw_bcid = word >> 11 & 0x1FF;  // Bits 11-19
        raw_time1 = word >> 6 & 0x1F;   // Bits 6-10
        raw_time2 = word >> 1 & 0x1F;   // Bits 1-5
    }
}

void correctBCWrapAround(int raw_bcid, int raw_bcout, int BC0, int& bcid, int& bcout) {
    bcid = raw_bcid - BC0;
    if (bcid < -128) bcid += 256;
    if (bcid > 128) bcid -= 256;

    bcout = raw_bcout % 256 - BC0;
    if (bcout < -128) bcout += 256;
    if (bcout > 128) bcout -= 256;
}

void storeRawHit(int clk, int word, int channel, int rise, int raw_bcid, int bcid, 
                 int raw_time1, int raw_time2, int raw_bcout, int bcout,
                 std::vector<int>& hit_clk, std::vector<int>& hit_channel, 
                 std::vector<int>& hit_raw_bcid, std::vector<int>& hit_bcid,
                 std::vector<int>& hit_time1, std::vector<int>& hit_time2, 
                 std::vector<int>& hit_rise, std::vector<int>& hit_raw_bcout, 
                 std::vector<int>& hit_bcout) {
    hit_clk.push_back(clk);
    hit_channel.push_back(channel);
    hit_raw_bcid.push_back(raw_bcid);
    hit_bcid.push_back(bcid);
    hit_time1.push_back(raw_time1);
    hit_time2.push_back(raw_time2);
    hit_rise.push_back(rise);
    hit_raw_bcout.push_back(raw_bcout);
    hit_bcout.push_back(bcout);
}

// ============================================================================================================
// Event processing functions

int extractTriggerTime(size_t n_hits, int trig_channel, 
                        const std::vector<int>& hit_channel, 
                        const std::vector<int>& hit_rise,
                        const std::vector<int>& hit_time2, 
                        const std::vector<int>& hit_bcid) {
    int trigger_time = -1;
    
    for (size_t idx_hit = 0; idx_hit < n_hits; idx_hit++) {
        if (hit_channel.at(idx_hit) == trig_channel && hit_rise.at(idx_hit) == 1 && trigger_time == -1) {
            if (hit_time2.at(idx_hit) != 0) {
                trigger_time = hit_bcid.at(idx_hit) * 30 + hit_time2.at(idx_hit) - 1;
            } else {
                std::cout << "WARNING: Trigger time empty" << std::endl;
            }
            break;
        }
    }
    
    return trigger_time;
}

void processHitsInEvent(size_t idx_hit, int n_hits, int trigger_time,
                        int dt_max, int dt_min, int trig_channel,
                        const std::vector<int>& hit_channel,
                        const std::vector<int>& hit_rise,
                        const std::vector<int>& hit_time1,
                        const std::vector<int>& hit_time2,
                        const std::vector<int>& hit_bcid,
                        std::vector<int>& proc_layer,
                        std::vector<int>& proc_strip,
                        std::vector<int>& proc_time1,
                        std::vector<int>& proc_time2,
                        std::vector<int>& proc_dt_time1_time2,
                        std::vector<int>& proc_dt_time1_trigger,
                        std::vector<int>& proc_dt_time2_trigger,
                        bool* eta1_hit_flags,
                        bool* eta2_hit_flags) {
    
    // Compute geometry
    int layer = (hit_channel.at(idx_hit) % 24) / 8;
    int column = hit_channel.at(idx_hit) / 24;
    int strip = 8 * column + hit_channel.at(idx_hit) % 8;
    
    // Compute time values
    int time1 = hit_time1.at(idx_hit) != 0 ? hit_bcid.at(idx_hit) * 30 + hit_time1.at(idx_hit) - 1 : -1;
    int time2 = hit_time2.at(idx_hit) != 0 ? hit_bcid.at(idx_hit) * 30 + hit_time2.at(idx_hit) - 1 : -1;
    
    // Efficiency counting for rising hits on non-trigger channels
    if (hit_rise.at(idx_hit) == 1 && hit_channel.at(idx_hit) != trig_channel) {
        if (time1 != -1) {
            int dt_time1_trigger = time1 - trigger_time;
            proc_dt_time1_trigger.push_back(dt_time1_trigger);
            
            if (dt_time1_trigger < dt_max && dt_time1_trigger > dt_min) {
                eta1_hit_flags[layer] = true;
            }
        }
        
        if (time2 != -1) {
            int dt_time2_trigger = time2 - trigger_time;
            proc_dt_time2_trigger.push_back(dt_time2_trigger);
            
            if (dt_time2_trigger < dt_max && dt_time2_trigger > dt_min) {
                eta2_hit_flags[layer] = true;
            }
        }
        
        if (time1 != -1 && time2 != -1) {
            int dt_time1_time2 = time1 - time2;
            proc_dt_time1_time2.push_back(dt_time1_time2);
        }
    }
    
    // Store processed hit data
    if (time1 != -1 || time2 != -1) {
        proc_layer.push_back(layer);
        proc_strip.push_back(strip);
    }
    if (time1 != -1) proc_time1.push_back(time1);
    if (time2 != -1) proc_time2.push_back(time2);
}

void clearEventVectors(std::vector<int>& hit_clk, std::vector<int>& hit_channel,
                       std::vector<int>& hit_raw_bcid, std::vector<int>& hit_bcid,
                       std::vector<int>& hit_time1, std::vector<int>& hit_time2,
                       std::vector<int>& hit_rise, std::vector<int>& hit_raw_bcout,
                       std::vector<int>& hit_bcout, std::vector<int>& proc_layer,
                       std::vector<int>& proc_strip, std::vector<int>& proc_time1,
                       std::vector<int>& proc_time2, std::vector<int>& proc_dt_time1_time2,
                       std::vector<int>& proc_trigger_time,
                       std::vector<int>& proc_dt_time1_trigger, std::vector<int>& proc_dt_time2_trigger,
                       std::vector<int>& proc_tot1, std::vector<int>& proc_tot2,
                       std::vector<int>& cluster_size1, std::vector<int>& cluster_size2,
                       std::vector<int>& cluster_tot1, std::vector<int>& cluster_tot2) {
    hit_clk.clear();
    hit_channel.clear();
    hit_raw_bcid.clear();
    hit_bcid.clear();
    hit_time1.clear();
    hit_time2.clear();
    hit_rise.clear();
    hit_raw_bcout.clear();
    hit_bcout.clear();
    
    proc_layer.clear();
    proc_strip.clear();
    proc_time1.clear();
    proc_time2.clear();
    proc_dt_time1_time2.clear();
    proc_trigger_time.clear();
    proc_dt_time1_trigger.clear();
    proc_dt_time2_trigger.clear();
    proc_tot1.clear();
    proc_tot2.clear();
    
    cluster_size1.clear();
    cluster_size2.clear();
    cluster_tot1.clear();
    cluster_tot2.clear();
}

void processEvent(int& n_event, int n_hits,
                  std::vector<int>& hit_clk, std::vector<int>& hit_channel, 
                  std::vector<int>& hit_raw_bcid, std::vector<int>& hit_bcid,
                  std::vector<int>& hit_time1, std::vector<int>& hit_time2, 
                  std::vector<int>& hit_rise, std::vector<int>& hit_raw_bcout, 
                  std::vector<int>& hit_bcout,
                  std::vector<int>& proc_layer, std::vector<int>& proc_strip,
                  std::vector<int>& proc_time1, std::vector<int>& proc_time2,
                  std::vector<int>& proc_dt_time1_time2,
                  std::vector<int>& proc_trigger_time,
                  std::vector<int>& proc_dt_time1_trigger, std::vector<int>& proc_dt_time2_trigger,
                  std::vector<int>& proc_tot1, std::vector<int>& proc_tot2,
                  std::vector<int>& cluster_size1, std::vector<int>& cluster_size2,
                  std::vector<int>& cluster_tot1, std::vector<int>& cluster_tot2,
                  int dt_max, int dt_min, int trig_channel,
                  EfficiencyCounters& counters) {
    
    n_event++;
    std::cout << "\nProcessing event No. " << n_event << std::endl;
    
    // Verify consistency
    if (n_hits != (int)hit_channel.size()) {
        std::cerr << "ERROR: n_hits mismatch!" << std::endl;
        return;
    }
    
    // Extract trigger time
    int trigger_time = extractTriggerTime(n_hits, trig_channel, hit_channel, hit_rise, hit_time2, hit_bcid);
    
    // Initialize local efficiency tracking for this event
    bool eta1_hit_flags[3] = {false, false, false};
    bool eta2_hit_flags[3] = {false, false, false};
    
    // Main hits processing loop
    for (size_t idx_hit = 0; idx_hit < (size_t)n_hits; idx_hit++) {
        processHitsInEvent(idx_hit, n_hits, trigger_time, dt_max, dt_min, trig_channel,
                          hit_channel, hit_rise, hit_time1, hit_time2, hit_bcid,
                          proc_layer, proc_strip, proc_time1, proc_time2,
                          proc_dt_time1_time2, proc_dt_time1_trigger, proc_dt_time2_trigger,
                          eta1_hit_flags, eta2_hit_flags);
    }
    
    // Store trigger time for event
    proc_trigger_time.push_back(trigger_time);
    
    // Update efficiency counters based on hit flags
    if (trigger_time != -1) {
        counters.triggered_events_external++;
        
        for (int layer = 0; layer < 3; layer++) {
            if (eta1_hit_flags[layer]) counters.eta1_efficiency_counter[layer]++;
            if (eta2_hit_flags[layer]) counters.eta2_efficiency_counter[layer]++;
            if (eta1_hit_flags[layer] || eta2_hit_flags[layer]) counters.eta_or_efficiency_counter[layer]++;
            if (eta1_hit_flags[layer] && eta2_hit_flags[layer]) counters.eta_and_efficiency_counter[layer]++;
        }
        
        // RPC trigger logic
        if ((eta1_hit_flags[1] && eta1_hit_flags[2]) || (eta2_hit_flags[1] && eta2_hit_flags[2])) {
            counters.triggered_events_rpc[0]++;
            if (eta1_hit_flags[0]) counters.eta1_efficiency_counter_rpc[0]++;
            if (eta2_hit_flags[0]) counters.eta2_efficiency_counter_rpc[0]++;
            if (eta1_hit_flags[0] || eta2_hit_flags[0]) counters.eta_or_efficiency_counter_rpc[0]++;
            if (eta1_hit_flags[0] && eta2_hit_flags[0]) counters.eta_and_efficiency_counter_rpc[0]++;
        }
        
        if ((eta1_hit_flags[0] && eta1_hit_flags[2]) || (eta2_hit_flags[0] && eta2_hit_flags[2])) {
            counters.triggered_events_rpc[1]++;
            if (eta1_hit_flags[1]) counters.eta1_efficiency_counter_rpc[1]++;
            if (eta2_hit_flags[1]) counters.eta2_efficiency_counter_rpc[1]++;
            if (eta1_hit_flags[1] || eta2_hit_flags[1]) counters.eta_or_efficiency_counter_rpc[1]++;
            if (eta1_hit_flags[1] && eta2_hit_flags[1]) counters.eta_and_efficiency_counter_rpc[1]++;
        }
        
        if ((eta1_hit_flags[0] && eta1_hit_flags[1]) || (eta2_hit_flags[0] && eta2_hit_flags[1])) {
            counters.triggered_events_rpc[2]++;
            if (eta1_hit_flags[2]) counters.eta1_efficiency_counter_rpc[2]++;
            if (eta2_hit_flags[2]) counters.eta2_efficiency_counter_rpc[2]++;
            if (eta1_hit_flags[2] || eta2_hit_flags[2]) counters.eta_or_efficiency_counter_rpc[2]++;
            if (eta1_hit_flags[2] && eta2_hit_flags[2]) counters.eta_and_efficiency_counter_rpc[2]++;
        }
    }
}

// ============================================================================================================
// Efficiency calculation functions

void calculateEfficiencies(const EfficiencyCounters& counters,
                          double* eta1_eff_ext, double* eta2_eff_ext, 
                          double* eta_or_eff_ext, double* eta_and_eff_ext,
                          double* eta1_eff_ext_err, double* eta2_eff_ext_err, 
                          double* eta_or_eff_ext_err, double* eta_and_eff_ext_err,
                          double* eta1_eff_rpc, double* eta2_eff_rpc, 
                          double* eta_or_eff_rpc, double* eta_and_eff_rpc,
                          double* eta1_eff_rpc_err, double* eta2_eff_rpc_err, 
                          double* eta_or_eff_rpc_err, double* eta_and_eff_rpc_err) {
    
    // External trigger only efficiency
    if (counters.triggered_events_external > 0) {
        for (int layer = 0; layer < 3; layer++) {
            eta1_eff_ext[layer] = static_cast<double>(counters.eta1_efficiency_counter[layer]) / counters.triggered_events_external;
            eta2_eff_ext[layer] = static_cast<double>(counters.eta2_efficiency_counter[layer]) / counters.triggered_events_external;
            eta_or_eff_ext[layer] = static_cast<double>(counters.eta_or_efficiency_counter[layer]) / counters.triggered_events_external;
            eta_and_eff_ext[layer] = static_cast<double>(counters.eta_and_efficiency_counter[layer]) / counters.triggered_events_external;
            
            eta1_eff_ext_err[layer] = std::sqrt(eta1_eff_ext[layer] * (1.0 - eta1_eff_ext[layer]) / counters.triggered_events_external);
            eta2_eff_ext_err[layer] = std::sqrt(eta2_eff_ext[layer] * (1.0 - eta2_eff_ext[layer]) / counters.triggered_events_external);
            eta_or_eff_ext_err[layer] = std::sqrt(eta_or_eff_ext[layer] * (1.0 - eta_or_eff_ext[layer]) / counters.triggered_events_external);
            eta_and_eff_ext_err[layer] = std::sqrt(eta_and_eff_ext[layer] * (1.0 - eta_and_eff_ext[layer]) / counters.triggered_events_external);
        }
    }
    
    // RPC trigger efficiency
    for (int layer = 0; layer < 3; layer++) {
        if (counters.triggered_events_rpc[layer] > 0) {
            eta1_eff_rpc[layer] = static_cast<double>(counters.eta1_efficiency_counter_rpc[layer]) / counters.triggered_events_rpc[layer];
            eta2_eff_rpc[layer] = static_cast<double>(counters.eta2_efficiency_counter_rpc[layer]) / counters.triggered_events_rpc[layer];
            eta_or_eff_rpc[layer] = static_cast<double>(counters.eta_or_efficiency_counter_rpc[layer]) / counters.triggered_events_rpc[layer];
            eta_and_eff_rpc[layer] = static_cast<double>(counters.eta_and_efficiency_counter_rpc[layer]) / counters.triggered_events_rpc[layer];
            
            eta1_eff_rpc_err[layer] = std::sqrt(eta1_eff_rpc[layer] * (1.0 - eta1_eff_rpc[layer]) / counters.triggered_events_rpc[layer]);
            eta2_eff_rpc_err[layer] = std::sqrt(eta2_eff_rpc[layer] * (1.0 - eta2_eff_rpc[layer]) / counters.triggered_events_rpc[layer]);
            eta_or_eff_rpc_err[layer] = std::sqrt(eta_or_eff_rpc[layer] * (1.0 - eta_or_eff_rpc[layer]) / counters.triggered_events_rpc[layer]);
            eta_and_eff_rpc_err[layer] = std::sqrt(eta_and_eff_rpc[layer] * (1.0 - eta_and_eff_rpc[layer]) / counters.triggered_events_rpc[layer]);
        }
    }
}

void createEfficiencyHistograms(TFile* output_file,
                               const double* eta1_eff_ext, const double* eta2_eff_ext,
                               const double* eta_or_eff_ext, const double* eta_and_eff_ext,
                               const double* eta1_eff_ext_err, const double* eta2_eff_ext_err,
                               const double* eta_or_eff_ext_err, const double* eta_and_eff_ext_err,
                               const double* eta1_eff_rpc, const double* eta2_eff_rpc,
                               const double* eta_or_eff_rpc, const double* eta_and_eff_rpc,
                               const double* eta1_eff_rpc_err, const double* eta2_eff_rpc_err,
                               const double* eta_or_eff_rpc_err, const double* eta_and_eff_rpc_err) {
    
    // External trigger histograms
    output_file->mkdir("efficiencies_histograms/external_trigger");
    output_file->cd("efficiencies_histograms/external_trigger");
    
    for (int layer = 0; layer < 3; layer++) {
        TH1F* hist = new TH1F(Form("eff_external_trigger%d", layer), 
                             Form("Layer %d Efficiency (External Trigger)", layer), 
                             4, 0.5, 4.5);
        hist->GetXaxis()->SetBinLabel(1, "eta1");
        hist->GetXaxis()->SetBinLabel(2, "eta2");
        hist->GetXaxis()->SetBinLabel(3, "OR");
        hist->GetXaxis()->SetBinLabel(4, "AND");
        
        hist->SetBinContent(1, eta1_eff_ext[layer]);
        hist->SetBinError(1, eta1_eff_ext_err[layer]);
        hist->SetBinContent(2, eta2_eff_ext[layer]);
        hist->SetBinError(2, eta2_eff_ext_err[layer]);
        hist->SetBinContent(3, eta_or_eff_ext[layer]);
        hist->SetBinError(3, eta_or_eff_ext_err[layer]);
        hist->SetBinContent(4, eta_and_eff_ext[layer]);
        hist->SetBinError(4, eta_and_eff_ext_err[layer]);
    }
    
    // RPC trigger histograms
    output_file->cd();
    output_file->mkdir("efficiencies_histograms/external_plus_rpc_trigger");
    output_file->cd("efficiencies_histograms/external_plus_rpc_trigger");
    
    for (int layer = 0; layer < 3; layer++) {
        TH1F* hist = new TH1F(Form("eff_rpc_layer%d", layer), 
                             Form("Layer %d Efficiency (External and RPC Trigger)", layer), 
                             4, 0.5, 4.5);
        hist->GetXaxis()->SetBinLabel(1, "eta1");
        hist->GetXaxis()->SetBinLabel(2, "eta2");
        hist->GetXaxis()->SetBinLabel(3, "OR");
        hist->GetXaxis()->SetBinLabel(4, "AND");
        
        hist->SetBinContent(1, eta1_eff_rpc[layer]);
        hist->SetBinError(1, eta1_eff_rpc_err[layer]);
        hist->SetBinContent(2, eta2_eff_rpc[layer]);
        hist->SetBinError(2, eta2_eff_rpc_err[layer]);
        hist->SetBinContent(3, eta_or_eff_rpc[layer]);
        hist->SetBinError(3, eta_or_eff_rpc_err[layer]);
        hist->SetBinContent(4, eta_and_eff_rpc[layer]);
        hist->SetBinError(4, eta_and_eff_rpc_err[layer]);
    }
}
