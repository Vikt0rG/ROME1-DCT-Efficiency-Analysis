#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <cstring>

#include "process_raw_hits.hpp"

// Parse command line arguments and return dt_max, dt_min, self_trigger
void parseArguments(int argc, char** argv, int& dt_max, int& dt_min, bool& self_trigger) {
    dt_max = -100;           // Default
    dt_min = -180;           // Default
    self_trigger = false;    // Default
    
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--dt-max") == 0 && i + 1 < argc) {
            dt_max = std::stoi(argv[++i]);
            std::cout << "Set dt_max = " << dt_max << std::endl;
        } else if (std::strcmp(argv[i], "--dt-min") == 0 && i + 1 < argc) {
            dt_min = std::stoi(argv[++i]);
            std::cout << "Set dt_min = " << dt_min << std::endl;
        } else if (std::strcmp(argv[i], "--self") == 0) {
            self_trigger = true;
            std::cout << "Using self-trigger mode" << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    std::cout << "Starting DCT raw hits processing..." << std::endl;
    
    // Parse command line arguments
    int dt_max;                             // Maximum time window for efficiency calculation
    int dt_min;                             // Minimum time window for efficiency calculation
    bool self_trigger;                      // Flag for self trigger condition
    parseArguments(argc, argv, dt_max, dt_min, self_trigger);
    const int trig_channel = 143;           // Trigger signal channel
    const int empty_word = 0x5555555;       // Empty word marker

    // Open input file
    std::ifstream tmp_file("tmp.strip");
    if (!tmp_file.is_open()) {
        std::cerr << "ERROR: Cannot open tmp.strip in current directory!" << std::endl;
        return 1;
    }

    // Setup output file and trees
    TFile* output_file = nullptr;
    TTree* input_data_tree = nullptr;
    TTree* processed_data_tree = nullptr;
    TTree* clusterization_tree = nullptr;
    
    setupOutputFile(output_file, input_data_tree, processed_data_tree, clusterization_tree);
    if (!output_file) {
        std::cerr << "ERROR: Output file setup failed!" << std::endl;
        return 1;
    }

    // Setup branch addresses
    std::vector<int> hit_clk, hit_channel, hit_raw_bcid, hit_bcid;
    std::vector<int> hit_time1, hit_time2, hit_rise, hit_raw_bcout, hit_bcout;
    std::vector<int> proc_layer, proc_strip, proc_time1, proc_time2;
    std::vector<int> proc_dt_time1_time2, proc_trigger_time;
    std::vector<int> proc_dt_time1_trigger, proc_dt_time2_trigger;
    std::vector<int> proc_tot1, proc_tot2;
    std::vector<int> cluster_size1, cluster_size2, cluster_tot1, cluster_tot2;

    int n_event = 0;
    int n_hits = 0;

    // Setup branches for raw data tree
    input_data_tree->Branch("hit_clk", &hit_clk);
    input_data_tree->Branch("hit_channel", &hit_channel);
    input_data_tree->Branch("hit_raw_bcid", &hit_raw_bcid);
    input_data_tree->Branch("hit_bcid", &hit_bcid);
    input_data_tree->Branch("hit_time1", &hit_time1);
    input_data_tree->Branch("hit_time2", &hit_time2);
    input_data_tree->Branch("hit_rise", &hit_rise);
    input_data_tree->Branch("hit_raw_bcout", &hit_raw_bcout);
    input_data_tree->Branch("hit_bcout", &hit_bcout);

    // Setup branches for processed data tree
    processed_data_tree->Branch("n_event", &n_event);
    processed_data_tree->Branch("n_hits", &n_hits);
    processed_data_tree->Branch("proc_layer", &proc_layer);
    processed_data_tree->Branch("proc_strip", &proc_strip);
    processed_data_tree->Branch("proc_time1", &proc_time1);
    processed_data_tree->Branch("proc_time2", &proc_time2);
    processed_data_tree->Branch("proc_dt_time1_time2", &proc_dt_time1_time2);
    processed_data_tree->Branch("trigger_time", &proc_trigger_time);
    processed_data_tree->Branch("proc_dt_time1_trigger", &proc_dt_time1_trigger);
    processed_data_tree->Branch("proc_dt_time2_trigger", &proc_dt_time2_trigger);
    processed_data_tree->Branch("proc_tot1", &proc_tot1);
    processed_data_tree->Branch("proc_tot2", &proc_tot2);

    // Setup branches for clusterization tree
    clusterization_tree->Branch("cluster_size1", &cluster_size1);
    clusterization_tree->Branch("cluster_size2", &cluster_size2);
    clusterization_tree->Branch("cluster_tot1", &cluster_tot1);
    clusterization_tree->Branch("cluster_tot2", &cluster_tot2);

    // Initialize efficiency counters
    EfficiencyCounters counters;
    initializeEfficiencyCounters(counters);

    // ============================================================================================================
    // Main processing loop
    int clk, word, raw_bcout;
    int BC0 = 0;

    while (tmp_file >> std::dec >> clk >> std::hex >> word >> raw_bcout) {
        
        // Skip empty words
        if (word == empty_word) continue;

        // Extract word data
        int channel, rise, raw_bcid, raw_time1, raw_time2;
        extractWordData(word, channel, rise, raw_bcid, raw_time1, raw_time2);

        // Initialize BC0 for new event
        if (n_hits == 0) BC0 = raw_bcid % 256;

        // Apply BC wrap-around correction
        int bcid, bcout;
        correctBCWrapAround(raw_bcid, raw_bcout, BC0, bcid, bcout);

        // Store raw hit
        storeRawHit(clk, word, channel, rise, raw_bcid, bcid, raw_time1, raw_time2, raw_bcout, bcout,
                   hit_clk, hit_channel, hit_raw_bcid, hit_bcid, hit_time1, hit_time2, 
                   hit_rise, hit_raw_bcout, hit_bcout);
        n_hits++;

        // WIP: For self-triggered data: Currently just skip events with too many hits
        if (self_trigger && n_hits > 20) {
            clearEventVectors(hit_clk, hit_channel, hit_raw_bcid, hit_bcid, hit_time1, hit_time2, 
                            hit_rise, hit_raw_bcout, hit_bcout, proc_layer, proc_strip, 
                            proc_time1, proc_time2, proc_dt_time1_time2, proc_trigger_time,
                            proc_dt_time1_trigger, proc_dt_time2_trigger, proc_tot1, proc_tot2,
                            cluster_size1, cluster_size2, cluster_tot1, cluster_tot2);
            n_hits = 0;
            continue;
        }

        // Process event at end of readout window (clk == 127)
        if (clk == 127) {
            processEvent(n_event, n_hits, hit_clk, hit_channel, hit_raw_bcid, hit_bcid,
                        hit_time1, hit_time2, hit_rise, hit_raw_bcout, hit_bcout,
                        proc_layer, proc_strip, proc_time1, proc_time2, proc_dt_time1_time2,
                        proc_trigger_time, proc_dt_time1_trigger, proc_dt_time2_trigger,
                        proc_tot1, proc_tot2, cluster_size1, cluster_size2, cluster_tot1, cluster_tot2,
                        dt_max, dt_min, trig_channel, counters);

            // Fill trees
            input_data_tree->Fill();
            processed_data_tree->Fill();
            clusterization_tree->Fill();

            // Clear vectors for next event
            clearEventVectors(hit_clk, hit_channel, hit_raw_bcid, hit_bcid, hit_time1, hit_time2, 
                            hit_rise, hit_raw_bcout, hit_bcout, proc_layer, proc_strip, 
                            proc_time1, proc_time2, proc_dt_time1_time2, proc_trigger_time,
                            proc_dt_time1_trigger, proc_dt_time2_trigger, proc_tot1, proc_tot2,
                            cluster_size1, cluster_size2, cluster_tot1, cluster_tot2);
            n_hits = 0;
        }
    }

    tmp_file.close();

    // ============================================================================================================
    // Calculate efficiencies and create histograms

    double eta1_eff_ext[3] = {0.0, 0.0, 0.0};
    double eta2_eff_ext[3] = {0.0, 0.0, 0.0};
    double eta_or_eff_ext[3] = {0.0, 0.0, 0.0};
    double eta_and_eff_ext[3] = {0.0, 0.0, 0.0};
    
    double eta1_eff_ext_err[3] = {0.0, 0.0, 0.0};
    double eta2_eff_ext_err[3] = {0.0, 0.0, 0.0};
    double eta_or_eff_ext_err[3] = {0.0, 0.0, 0.0};
    double eta_and_eff_ext_err[3] = {0.0, 0.0, 0.0};

    double eta1_eff_rpc[3] = {0.0, 0.0, 0.0};
    double eta2_eff_rpc[3] = {0.0, 0.0, 0.0};
    double eta_or_eff_rpc[3] = {0.0, 0.0, 0.0};
    double eta_and_eff_rpc[3] = {0.0, 0.0, 0.0};

    double eta1_eff_rpc_err[3] = {0.0, 0.0, 0.0};
    double eta2_eff_rpc_err[3] = {0.0, 0.0, 0.0};
    double eta_or_eff_rpc_err[3] = {0.0, 0.0, 0.0};
    double eta_and_eff_rpc_err[3] = {0.0, 0.0, 0.0};

    calculateEfficiencies(counters, eta1_eff_ext, eta2_eff_ext, eta_or_eff_ext, eta_and_eff_ext,
                         eta1_eff_ext_err, eta2_eff_ext_err, eta_or_eff_ext_err, eta_and_eff_ext_err,
                         eta1_eff_rpc, eta2_eff_rpc, eta_or_eff_rpc, eta_and_eff_rpc,
                         eta1_eff_rpc_err, eta2_eff_rpc_err, eta_or_eff_rpc_err, eta_and_eff_rpc_err);

    createEfficiencyHistograms(output_file,
                              eta1_eff_ext, eta2_eff_ext, eta_or_eff_ext, eta_and_eff_ext,
                              eta1_eff_ext_err, eta2_eff_ext_err, eta_or_eff_ext_err, eta_and_eff_ext_err,
                              eta1_eff_rpc, eta2_eff_rpc, eta_or_eff_rpc, eta_and_eff_rpc,
                              eta1_eff_rpc_err, eta2_eff_rpc_err, eta_or_eff_rpc_err, eta_and_eff_rpc_err);

    // Write output and cleanup
    output_file->Write();
    output_file->Close();

    std::cout << "Processing complete. Output written to out.root" << std::endl;
    return 0;
}