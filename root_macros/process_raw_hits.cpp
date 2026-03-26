#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#include <fstream>
#include <iostream>

#include <vector>

// #include "utils.hpp"


void process_raw_hits() {

    // TODO: Just provide these two as arguments from shell skript and also store dt_hit_time*_trigger to the output for user to easily readjust time window
    int dt_max = -100;
    int dt_min = -180;

    // TODO: Add self-trigger option (if True, skip events that have more than 20 hits)
    bool self_trigger = true;

    std::ifstream tmp_file("tmp.strip");
    if (!tmp_file.is_open()) {
        std::cerr << "ERROR: Cannot open tmp.strip in current directory!" << std::endl;
        return;
    }

    // Create output ROOT file and define trees
    TFile* output_file = TFile::Open("out.root", "RECREATE");
    if (!output_file) return;
    TTree* input_data_tree = new TTree("input_data", "Input DCT data");
    TTree* processed_data_tree = new TTree("processed_data", "Processed DCT data");
    TTree* clusterization_tree = new TTree("clusterization_data", "Data after clusterization and corrections");


    // ============================================================================================================
    // Global constants
    const int trig_channel = 143;                  // Trigger signal from scintillator is treated as a DCT hit and is sent to channel 143


    // TODO: Update variable initialization: Name conventions for different types, add documentation etc etc.

    // Global variables and set initial counters
    int n_event = 0;
    int n_hits = 0;
    int n_processed_hits = 0;
    int n_clusters_eta1 = 0;
    int n_clusters_eta2 = 0;
    int cluster_size_eta1 = 0;
    int cluster_size_eta2 = 0;
    int n_tracks_eta1 = 0;
    int n_tracks_eta2 = 0;
    int track_length_eta1 = 0;
    int track_length_eta2 = 0;
    int trigger_time = -1;                         // Initialize trigger time to -1 (invalid value) to identify events without a trigger hit

    // Vectors for the DCT output: (clk, word, raw_bcout)
    std::vector<int> hit_clk;                      // Time of hit in readout window (tick = 3.125 ns)

    // Word structure (taken from DCT 'documentation' at https://gitlab.cern.ch/groups/atlas-tdaq-muon-barrel-trigger/-/wikis/home/DCT/%7BBI-DCT-test-station-at-BB5%7D):
    std::vector<int> hit_channel;                  // 8-bit strip ID
    std::vector<int> hit_raw_bcid;                 // 8-bit (or 9-bit for falling edge) BC ID
    std::vector<int> hit_bcid;                     // 8-bit (or 9-bit for falling edge) BC ID after BC wrap-around correction
    std::vector<int> hit_time1;                    // 5-bit time of η1
    std::vector<int> hit_time2;                    // 5-bit time of η2
    std::vector<int> hit_rise;                     // 1-bit rise or fall 'status' (1 = rise, 0 = fall)

    std::vector<int> hit_raw_bcout;                // Number of BC at which word exits the DCT
    std::vector<int> hit_bcout;                    // Number of BC at which word exits the DCT after BC wrap-around correction

    // Vectors for processed data
    std::vector<int> proc_layer;                   // Detector layer (0, 1 or 2) of the hit (stored only if either time1 OR time2 is present)
    std::vector<int> proc_strip;                   // Strip number (0-47) of the hit (stored only if either time1 OR time2 is present)
    std::vector<int> proc_time1, proc_time2;       // Time of hit in η1 and η2 (stored only if hit_time* != 0), converted from raw time and BCID information, units of ticks (0.833 ns)
    std::vector<int> proc_dt_time1_time2;          // dt between time1 and time2 for each hit (only when both times measured), units of ticks (0.833 ns)
    std::vector<int> proc_trigger_time;            // Calculated trigger time for each event (same for all hits in the same event), units of ticks (0.833 ns)
    std::vector<int> proc_dt_time1_trigger,
                     proc_dt_time2_trigger;        // dt between hit time and trigger time for efficiency calculation, units of ticks (0.833 ns)
    std::vector<int> proc_tot1, proc_tot2;         // dt of the rising hit and corresponding falling hit (only if present), units of ticks (0.833 ns)

    // Vectors for data after clusterization
    std::vector<int> v_cluster_size_eta1, v_cluster_size_eta2; // Number of consecutive strips firing within 18 ticks (15 ns) time window within the same event
    std::vector<int> v_cluster_tot1, v_cluster_tot2;   // ToT calculated only for cluster centers (first hits in clusters)

    // Vectors for track reconstruction data (WIP)
    std::vector<int> v_track_length_eta1, v_track_length_eta2; // Number of consecutive layers with hits in the same event


    // ============================================================================================================
    // Branch definitions for raw data tree
    input_data_tree->Branch("hit_clk", &hit_clk);
    input_data_tree->Branch("hit_channel", &hit_channel);
    input_data_tree->Branch("hit_raw_bcid", &hit_raw_bcid);
    input_data_tree->Branch("hit_bcid", &hit_bcid);
    input_data_tree->Branch("hit_time1", &hit_time1);
    input_data_tree->Branch("hit_time2", &hit_time2);
    input_data_tree->Branch("hit_rise", &hit_rise);
    input_data_tree->Branch("hit_raw_bcout", &hit_raw_bcout);
    input_data_tree->Branch("hit_bcout", &hit_bcout);

    // Branch definitions for processed data tree
    processed_data_tree->Branch("n_event", &n_event);
    processed_data_tree->Branch("n_hits", &n_hits);
    processed_data_tree->Branch("proc_layer", &proc_layer);
    processed_data_tree->Branch("proc_strip", &proc_strip);
    processed_data_tree->Branch("proc_time1", &proc_time1);
    processed_data_tree->Branch("proc_time2", &proc_time2);
    processed_data_tree->Branch("proc_dt_time1_time2",&proc_dt_time1_time2);
    processed_data_tree->Branch("trigger_time", &proc_trigger_time);
    processed_data_tree->Branch("proc_dt_time1_trigger",&proc_dt_time1_trigger);
    processed_data_tree->Branch("proc_dt_time2_trigger",&proc_dt_time2_trigger);
    processed_data_tree->Branch("proc_tot1",&proc_tot1);
    processed_data_tree->Branch("proc_tot2",&proc_tot2);

    // Branch definitions for clusterization tree
    clusterization_tree->Branch("cluster_size_eta1", &cluster_size_eta1);
    clusterization_tree->Branch("cluster_size_eta2", &cluster_size_eta2);
    clusterization_tree->Branch("v_cluster_tot1", &v_cluster_tot1);
    clusterization_tree->Branch("v_cluster_tot2", &v_cluster_tot2);


    // ============================================================================================================
    // Initialize (raw) variables for words and hits processing loop.
    int clk, word, raw_bcout;
    int raw_bcid, BC0, bcid, bcout, raw_time1, raw_time2;

    // Efficiency counters
    int triggered_events_external = 0;                  // Counter for number of triggered events with an external trigger
    int triggered_events_rpc[3] = {0, 0, 0};            // Counter for number of triggered events with an RPC as a trigger

    // External trigger only efficiency counter arrays
    int eta1_efficiency_counter[3] = {0, 0, 0};         // Efficiency counters for each η layer
    int eta2_efficiency_counter[3] = {0, 0, 0};
    int eta_or_efficiency_counter[3] = {0, 0, 0};       // Efficiency counters for OR of both sides in each η layer
    int eta_and_efficiency_counter[3] = {0, 0, 0};      // Efficiency counters for AND of both sides in each η layer

    // External trigger + RPC as a trigger efficiency counter arrays
    int eta1_efficiency_counter_rpc[3] = {0, 0, 0};     // Efficiency counters for each η layer
    int eta2_efficiency_counter_rpc[3] = {0, 0, 0};
    int eta_or_efficiency_counter_rpc[3] = {0, 0, 0};   // Efficiency counters for OR of both sides in each η layer
    int eta_and_efficiency_counter_rpc[3] = {0, 0, 0};  // Efficiency counters for AND of both sides in each η layer

    // Start main loop to process raw data/words.
    while (tmp_file >> std::dec >> clk >> std::hex >> word >> raw_bcout) {
    
        // Local counters and flags
        int eta1_hits_counter[3] = {0, 0, 0};           // Counters for number of hits in each η layer (for later analysis)
        int eta2_hits_counter[3] = {0, 0, 0};
        bool eta1_hit_flags[3] = {false, false, false}; // Flags to identify if at least one hit is present in each η layer (for efficiency calculation)
        bool eta2_hit_flags[3] = {false, false, false};

        // Extract information from the word
        int channel = (word >> 20) & 0xFF;    // Bits 20-27
        int rise = word & 0x01;               // Bit 0

        if (rise == 1) {                      // Rising edge
            raw_bcid = word >> 12 & 0xFF;       // Bits 12-19
            raw_time1 = word >> 7 & 0x1F;       // Bits 7-11
            raw_time2 = word >> 1 & 0x3F;       // Bits 1-5

        } else {                              // Falling edge
            raw_bcid = word >> 11 & 0x1FF;      // Bits 11-19
            raw_time1 = word >> 6 & 0x1F;       // Bits 6-10
            raw_time2 = word >> 1 & 0x1F;       // Bits 1-5
        }

        // Process raw event information (skip empty words with value 0x5555555)
        if (word != 0x5555555) {
            // For each new event, n_hits is reset back to 0 and a new BC0 is recomputed
            if (n_hits == 0) BC0 = raw_bcid % 256; // Recompute reference BCID (BC0) for each new event
            // std::cout << BC0 << " " << n_hits << std::endl; // Debug print to check if BC0 is computed correctly
            n_hits++;

            // Handle BC wrap-around
            bcid = raw_bcid - BC0;
            if (bcid < -128) bcid += 256;
            if (bcid > 128) bcid -= 256;

            bcout = raw_bcout % 256 - BC0;
            if (bcout < -128) bcout += 256;
            if (bcout > 128) bcout -= 256;

            /* INSERT: BC wrap-around
            Initial ranges: [0, 255] for rising edge and [0, 511] for falling edge
            Step 1: raw_bcid - BC0          // Range: [-255, 255]
            Step 2: +/- 128                 // Range: [-128, 128]
            */

            // Store hit information in vectors
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

        // Skip to the next event if self-trigger option is enabled and more than one hit is present in the current event
        if (self_trigger && n_hits > 20) {
            // Clear hit vectors for the next event
            hit_clk.clear();

            hit_channel.clear();
            hit_raw_bcid.clear();
            hit_bcid.clear();
            hit_time1.clear();
            hit_time2.clear();
            hit_rise.clear();

            hit_raw_bcout.clear();
            hit_bcout.clear();

            n_hits = 0; // Reset hit counter for the next event
            continue;   // Skip to the next event
        }

        // --------------------------------------------------------------------------------------------------------
        // Process end of event (when clk reaches 127, which corresponds to the end of the 128-tick readout window)
        if (clk == 127) {
            n_event++;
            std::cout << "\n" << std::endl;
            std::cout << "==================================================================" << std::endl;
            std::cout << "Processing event No. " << n_event << std::endl;

            // Verify n_hits
            if (n_hits != hit_channel.size()) {
                std::cerr << "ERROR: n_hits " << n_hits << " does not match size of hit_channel vector " << hit_channel.size() << "!" << std::endl;
                return;
            }

            // Trigger time extraction loop (look for rising hit on trigger channel and extract trigger time)
            for (size_t idx_hit = 0; idx_hit < n_hits; idx_hit++) {
                if (hit_channel.at(idx_hit) == trig_channel && hit_rise.at(idx_hit) == 1 && trigger_time == -1) { 
                    // Trigger signal placed on channel 287 (143th channel on the eta2 side -> Check for hit/time2 != 0)
                    if (hit_time2.at(idx_hit) != 0) {
                        // NOTE: All time values are expressed in units of 0.833 ns ticks (30 ticks in one BC)
                        // A -1 subtraction is done to hit_time since hit_time starts from 1 (0 := no edge was measured)
                        trigger_time = hit_bcid.at(idx_hit) * 30 + hit_time2.at(idx_hit) - 1;
                    } else {
                        std::cout << "WARNING: Trigger time empty" << std::endl;
                    }
	            }
            }

            // ----------------------------------------------------------------------------------------------------
            // Hits loop and calculations section (e.g. efficiency counting, ToT calculation, etc.)

            // TODO: Geometry and time arrays for rising hits only
            std::vector<int> hit_layer(n_hits, -1);             // Detector layer (0, 1 or 2) of the hit (initialized to -1 for hits without time information)
            std::vector<int> hit_strip(n_hits, -1);             // Strip number (0-47) of the hit (initialized to -1 for hits without time information)
            std::vector<int> hit_time1_converted(n_hits, -1);   // Time of hit in η1 (if present, otherwise -1), converted from raw time and BCID information
            std::vector<int> hit_time2_converted(n_hits, -1);   // Time of hit in η2 (if present, otherwise -1), converted from raw time and BCID information

            // First hits loop to process geometry and time information of rising hits only (to be used in clusterization and track reconstruction)
            for (size_t idx_hit = 0; idx_hit < n_hits; idx_hit++) {

                if (hit_rise.at(idx_hit) == 0) continue;                // Skip falling hits

                // Detector geometry variables
                int layer = (hit_channel.at(idx_hit) % 24) / 8;         // Detector layer (0, 1 or 2)
                int column = hit_channel.at(idx_hit) / 24;              // Detector column (0, 1, 2, ... 5)
                int strip = 8 * column + hit_channel.at(idx_hit) % 8;   // Strip number (0-47)

                // Time variables to be used in efficiency counting and ToT calculation
                int time1 = hit_time1.at(idx_hit) != 0 ? hit_bcid.at(idx_hit) * 30 + hit_time1.at(idx_hit) - 1 : -1; // Time of hit in η1 (if present, otherwise -1)
                int time2 = hit_time2.at(idx_hit) != 0 ? hit_bcid.at(idx_hit) * 30 + hit_time2.at(idx_hit) - 1 : -1; // Time of hit in η2 (if present, otherwise -1)

                // Store geometry and time information in the hit_layer, hit_strip, hit_time*_converted arrays for later use in clusterization and track reconstruction
                hit_layer.at(idx_hit) = layer;
                hit_strip.at(idx_hit) = strip;
                hit_time1_converted.at(idx_hit) = time1;
                hit_time2_converted.at(idx_hit) = time2;

                // Push back 
                if (time1 != -1 || time2 != -1) proc_layer.push_back(layer);
                if (time1 != -1 || time2 != -1) proc_strip.push_back(strip);
                if (time1 != -1) proc_time1.push_back(time1);
                if (time2 != -1) proc_time2.push_back(time2);

            } // End of hits loop #1 processing geometry and time information of rising hits only

            // Clusterization logic variables
            std::vector<int> v_cluster_hit_assinged_eta1(n_hits, -1);    // -1 = not set, other values correspond to the index of the cluster the hit is assigned to
            std::vector<int> v_cluster_hit_assinged_eta2(n_hits, -1);
            std::vector<int> v_cluster_hit_is_center_eta1(n_hits, -1);   // -1 = not set, 0 = not cluster center, 1 = cluster center (first in time hit in cluster)
            std::vector<int> v_cluster_hit_is_center_eta2(n_hits, -1);

            // Vectors to store information of hits belonging to the same cluster
            std::vector<int> v_cluster_strip_eta1, v_cluster_strip_eta2; // Strip numbers of hits belonging to the same cluster (to check if they are consecutive)
            std::vector<int> v_cluster_time_eta1, v_cluster_time_eta2;   // Rising times of hits belonging to the same cluster (to check if they are within a certain time window from each other)

            // Track reconstruction logic arrays
            std::vector<int> track_hit_assinged_eta1(n_hits, -1);        // -1 = hit not assigned to a track, other values correspond to the index of the track the hit is assigned to
            std::vector<int> track_hit_assinged_eta2(n_hits, -1);

            // Vectors to store information of hits belonging to the same track
            std::vector<int> v_track_strip_eta1, v_track_strip_eta2;     // Strip numbers of hits belonging to the same track
            std::vector<int> v_track_time_eta1, v_track_time_eta2;       // Times of hits belonging to the same track

            // Arrays for ToT calculation after clusterization (only for cluster centers)
            std::vector<int> cluster_hit_used_tot_eta1(n_hits, -1);      // -1 = initial value or no information, 0 = hit has been used as falling hit, 1 = hit has been used as rising hit 
            std::vector<int> cluster_hit_used_tot_eta2(n_hits, -1);

            // TODO: Decide whether to keep ToT calculation before doing clusterization at all
            // Arrays for ToT calculation (before clusterization)
            std::vector<int> hit_used_tot_eta1(n_hits, -1);              // -1 = initial value or no information, 0 = hit has been used as falling hit, 1 = hit has been used as rising hit
            std::vector<int> hit_used_tot_eta2(n_hits, -1);

            // Second hits loop to process all hits and do clusterization, track reconstruction, efficiency counting, ToT calculation, etc.
            for (size_t idx_hit = 0; idx_hit < n_hits; idx_hit++) {

                // Detector geometry variables
                int layer = hit_layer.at(idx_hit);                      // Detector layer (0, 1 or 2)
                int strip = hit_strip.at(idx_hit);                      // Strip number (0-47)

                // Time variables to be used in efficiency counting and ToT calculation
                int time1 = hit_time1_converted.at(idx_hit);            // Time of hit in η1 (if present, otherwise -1)
                int time2 = hit_time2_converted.at(idx_hit);            // Time of hit in η2 (if present, otherwise -1)
                int tot1 = -1;                                          // Initial Time over Threshold (η1)
                int tot2 = -1;                                          // Initial Time over Threshold (η2)

                // Clusterization and track reconstruction logics as well as efficiency counting and ToT calculation for each hit in the event
                if (hit_rise.at(idx_hit) == 1 && hit_channel.at(idx_hit) != trig_channel) {

                    // --------------------------------------------------------------------------------------------
                    // Clusterization logic subsection

                    // For side η1
                    std::cout << "------------------------------------------------------------------" << std::endl;
                    std::cout << "Processing side η1: Hit: " << idx_hit << "; Layer " << layer << "; Strip: " << strip << "; Time " << time1 << std::endl;

                    // Initialize a new cluster containing only the current hit if the hit is not already in a cluster and time information is present
                    if (v_cluster_hit_assinged_eta1.at(idx_hit) == -1 && time1 != -1) {

                        v_cluster_strip_eta1.clear();                       // Clear cluster_strip* vector for the new cluster
                        v_cluster_time_eta1.clear();                        // Clear cluster_time* vector for the new cluster

                        v_cluster_strip_eta1.push_back(strip);              // Push back the strip number of the current hit for later checks of strip adjacency between cluster members
                        v_cluster_time_eta1.push_back(time1);               // Push back the time of the current hit for later time window checks between cluster members and potential cluster partner hits

                        v_cluster_hit_assinged_eta1.at(idx_hit) = n_clusters_eta1;    // Mark current hit as belonging to a cluster
                        v_cluster_hit_is_center_eta1.at(idx_hit) = 1;       // Define current hit as initial cluster center (to be updated later if another hit in the cluster has smaller rising time)
                        int idx_hit_cluster_center = idx_hit;             // Define current hit as initial cluster center (to be updated later if another hit in the cluster has smaller rising time)

                        n_clusters_eta1++;
                        cluster_size_eta1++;

                        // Iterate through the following hits of the same event and look for the potential cluster partners
                        for (size_t idx_next_hit = idx_hit + 1; idx_next_hit < n_hits; idx_next_hit++) {
                            if (v_cluster_hit_assinged_eta1.at(idx_next_hit) == -1 &&
                                hit_time1.at(idx_next_hit) != -1 &&
                                hit_rise.at(idx_next_hit) == 1 && 
                                hit_channel.at(idx_next_hit) != trig_channel) { // Only consider rising hits on non-trigger channels with valid time information as potential cluster partners

                                int next_layer = hit_layer.at(idx_next_hit);
                                if (next_layer != layer) continue;

                                // TODO: Decide whether to introduce any kind of metric in time and strip for clusterization and track reconstruction
                                // Only consider hits on the same or neighboring strips as potential cluster partners
                                bool strip_adjacent = false;
                                int next_strip = hit_strip.at(idx_next_hit);
                                for (size_t idx_cluster_hit = 0; idx_cluster_hit < cluster_size_eta1; idx_cluster_hit++) {
                                    if (abs(next_strip - v_cluster_strip_eta1.at(idx_cluster_hit)) <= 1) {
                                        // If at least one hit in the cluster is on the same or neighboring strip as the potential cluster partner hit, we can consider this potential cluster partner hit for further checks
                                        strip_adjacent = true;
                                        break;
                                    }
                                }
                                if (!strip_adjacent) continue; // If none of the hits in the cluster are on the same or neighboring strip as the potential cluster partner hit, skip this potential cluster partner hits
                                std::cout << "  Potential cluster partner: Hit " << idx_next_hit << "; Layer: " << next_layer << "; Strip: " << next_strip << "; Time " << hit_time1.at(idx_next_hit) << std::endl;

                                // Check if the potential cluster partner hit is within a certain time window (e.g. 15 ns = 18 ticks) from any of the hits already belonging to the cluster
                                int next_time1 = hit_time1_converted.at(idx_next_hit);
                                if (next_time1 == -1) continue;
                                for (size_t idx_cluster_hit = 0; idx_cluster_hit < cluster_size_eta1; idx_cluster_hit++) {
                                    if (abs(next_time1 - v_cluster_time_eta1.at(idx_cluster_hit)) < 18) {
                                        v_cluster_hit_assinged_eta1.at(idx_next_hit) = n_clusters_eta1;    // Mark the new cluster partner hit as belonging to the same cluster
                                        v_cluster_strip_eta1.push_back(next_strip); // Push back the strip number of the new cluster member hit
                                        v_cluster_time_eta1.push_back(next_time1);  // Push back the time of the new cluster member hit
                                        cluster_size_eta1++; // Increase cluster size

                                        // Check wether the new cluster partner hit has smaller rising time than the current cluster center hit, and if so, update the cluster center and its time
                                        if (next_time1 < time1) {
                                            idx_hit_cluster_center = idx_next_hit;
                                            v_cluster_hit_is_center_eta1.at(idx_hit_cluster_center) = 1; // Mark new cluster partner hit as cluster center
                                            v_cluster_hit_is_center_eta1.at(idx_hit) = 0; // Mark old cluster center hit as non-cluster center
                                            time1 = next_time1; // Update cluster center time to the smaller one between the current cluster center and the new cluster partner hit
                                        }
                                        break; // No need to check other hits in the v_cluster_time_eta1 vector, since we just need at least one hit in the cluster to be within the time window from the potential cluster partner hit
                                        // WIP: breaking here is still questionable since this way clusters can be very "long" in time
                                    }
                                }

                            }
                        } // End of loop looking for potential cluster partners

                        // If no more potential cluster partners are found for the current hit, store cluster size and reset cluster size variable for the next cluster
                        v_cluster_size_eta1.push_back(cluster_size_eta1);
                        n_clusters_eta1 = 0; // Reset cluster counter for the next cluster

                        // DEBUG: Print out parameters of the hits belonging to the cluster for the current hit
                        std::cout << "******************************************************************" << std::endl;
                        std::cout << "Resulting cluster: " << std::endl;
                        for (size_t idx_cluster_hit = 0; idx_cluster_hit < v_cluster_time_eta1.size(); idx_cluster_hit++) {
                            std::cout << idx_cluster_hit << ": strip = " << v_cluster_strip_eta1.at(idx_cluster_hit) << ", time = " << v_cluster_time_eta1.at(idx_cluster_hit) << std::endl;
                        }
                        std::cout << "\n";

                    }   // End of finding cluster center and its partners for side η1

                    // For side η2

                    // Initialize a new cluster containing only the current hit if the hit is not already in a cluster and time information is present
                    if (v_cluster_hit_assinged_eta2.at(idx_hit) == -1 && time2 != -1) {

                        v_cluster_strip_eta2.clear();                       // Clear cluster_strip* vector for the new cluster
                        v_cluster_time_eta2.clear();                        // Clear cluster_time* vector for the new cluster

                        v_cluster_strip_eta2.push_back(strip);              // Push back the strip number of the current hit for later checks of strip adjacency between cluster members
                        v_cluster_time_eta2.push_back(time2);               // Push back the time of the current hit for later time window checks between cluster members and potential cluster partner hits

                        v_cluster_hit_assinged_eta2.at(idx_hit) = n_clusters_eta2;    // Mark current hit as belonging to a cluster
                        v_cluster_hit_is_center_eta2.at(idx_hit) = 1;       // Define current hit as initial cluster center (to be updated later if another hit in the cluster has smaller rising time)
                        int idx_hit_cluster_center = idx_hit;             // Define current hit as initial cluster center (to be updated later if another hit in the cluster has smaller rising time)

                        n_clusters_eta2++;
                        cluster_size_eta2++;

                        // Iterate through the following hits of the same event and look for the potential cluster partners
                        for (size_t idx_next_hit = idx_hit + 1; idx_next_hit < n_hits; idx_next_hit++) {
                            if (v_cluster_hit_assinged_eta2.at(idx_next_hit) == -1 &&
                                hit_time2.at(idx_next_hit) != -1 &&
                                hit_rise.at(idx_next_hit) == 1 && 
                                hit_channel.at(idx_next_hit) != trig_channel) { // Only consider rising hits on non-trigger channels with valid time information as potential cluster partners

                                int next_layer = hit_layer.at(idx_next_hit);
                                if (next_layer != layer) continue;

                                // TODO: Decide whether to introduce any kind of metric in time and strip for clusterization and track reconstruction
                                // Only consider hits on the same or neighboring strips as potential cluster partners
                                bool strip_adjacent = false;
                                int next_strip = hit_strip.at(idx_next_hit);
                                for (size_t idx_cluster_hit = 0; idx_cluster_hit < cluster_size_eta2; idx_cluster_hit++) {
                                    if (abs(next_strip - v_cluster_strip_eta2.at(idx_cluster_hit)) <= 1) {
                                        // If at least one hit in the cluster is on the same or neighboring strip as the potential cluster partner hit, we can consider this potential cluster partner hit for further checks
                                        strip_adjacent = true;
                                        break;
                                    }
                                }
                                if (!strip_adjacent) continue; // If none of the hits in the cluster are on the same or neighboring strip as the potential cluster partner hit, skip this potential cluster partner hits

                                // Check if the potential cluster partner hit is within a certain time window (e.g. 15 ns = 18 ticks) from any of the hits already belonging to the cluster
                                int next_time2 = hit_time2_converted.at(idx_next_hit);
                                if (next_time2 == -1) continue;
                                for (size_t idx_cluster_hit = 0; idx_cluster_hit < cluster_size_eta2; idx_cluster_hit++) {
                                    if (abs(next_time2 - v_cluster_time_eta2.at(idx_cluster_hit)) < 18) {
                                        v_cluster_hit_assinged_eta2.at(idx_next_hit) = n_clusters_eta2;    // Mark the new cluster partner hit as belonging to the same cluster
                                        v_cluster_strip_eta2.push_back(next_strip); // Push back the strip number of the new cluster member hit
                                        v_cluster_time_eta2.push_back(next_time2);  // Push back the time of the new cluster member hit
                                        cluster_size_eta2++; // Increase cluster size

                                        // Check wether the new cluster partner hit has smaller rising time than the current cluster center hit, and if so, update the cluster center and its time
                                        if (next_time2 < time2) {
                                            idx_hit_cluster_center = idx_next_hit;
                                            v_cluster_hit_is_center_eta2.at(idx_hit_cluster_center) = 1; // Mark new cluster partner hit as cluster center
                                            v_cluster_hit_is_center_eta2.at(idx_hit) = 0; // Mark old cluster center hit as non-cluster center
                                            time2 = next_time2; // Update cluster center time to the smaller one between the current cluster center and the new cluster partner hit
                                        }
                                        break; // No need to check other hits in the v_cluster_time_eta2 vector, since we just need at least one hit in the cluster to be within the time window from the potential cluster partner hit
                                        // WIP: breaking here is still questionable since this way clusters can be very "long" in time
                                    }
                                }

                            }
                        } // End of loop looking for potential cluster partners

                        // If no more potential cluster partners are found for the current hit, store cluster size and reset cluster size variable for the next cluster
                        v_cluster_size_eta2.push_back(cluster_size_eta2);
                        n_clusters_eta2 = 0;

                    }   // End of finding cluster center and its partners for side η2

                    // --------------------------------------------------------------------------------------------
                    // WIP: Track reconstruction logic subsection

                    // For side η1

                    // Initialize a new track containing only the current hit if the hit is not already in a track and time information is present
                    if (track_hit_assinged_eta1.at(idx_hit) == -1 && time1 != -1) {

                        v_track_strip_eta1.clear();                       // Clear track_strip* vector for the new track
                        v_track_time_eta1.clear();                        // Clear track_time* vector for the new track

                        v_track_strip_eta1.push_back(strip);              // Push back the strip number of the current hit for later checks of strip adjacency between track members
                        v_track_time_eta1.push_back(time1);               // Push back the time of the current hit for later time window checks between track members and potential track partner hits

                        track_hit_assinged_eta1.at(idx_hit) = n_tracks_eta1;    // Mark current hit as belonging to a track

                        n_tracks_eta1++;
                        track_length_eta1++;

                        // Iterate through the following hits of the same event and look for the potential track partners
                        for (size_t idx_next_hit = idx_hit + 1; idx_next_hit < n_hits; idx_next_hit++) {
                            if (track_hit_assinged_eta1.at(idx_next_hit) == -1 &&
                                hit_time1.at(idx_next_hit) != -1 &&
                                hit_rise.at(idx_next_hit) == 1 && 
                                hit_channel.at(idx_next_hit) != trig_channel) { // Only consider rising hits on non-trigger channels with valid time information as potential track partners

                                int next_layer = hit_layer.at(idx_next_hit);
                                if (next_layer == layer) continue; // Only consider hits on different layers as potential track partners


                        // If no more potential track partners are found for the current hit, store track length and reset track length variable for the next track
                        v_track_length_eta1.push_back(track_length_eta1);
                        track_length_eta1 = 0;

                    } // End of track reconstruction logic for side η1

                    // --------------------------------------------------------------------------------------------
                    // Efficiency counters update subsection
                    // TODO: Move efficiency calculation after clusterization and track reconstruction
                    // TODO: Require hits from close-lying strips (add distance metric as an argument)
                    if (time1 != -1) {
                        int dt_time1_trigger = time1 - trigger_time;
                        proc_dt_time1_trigger.push_back(dt_time1_trigger);
                        hit_used_tot_eta1.at(idx_hit) = 1;      // Mark hit as used as rising hit for η1

                        if (dt_time1_trigger < dt_max &&
                            dt_time1_trigger > dt_min) {
                                eta1_hits_counter[layer]++;
                                if (eta1_hit_flags[layer] == false) {
                                    eta1_hit_flags[layer] = true;
                                }
                            }
                    }
                    if (time2 != -1) {
                        int dt_time2_trigger = time2 - trigger_time;
                        proc_dt_time2_trigger.push_back(dt_time2_trigger);
                        hit_used_tot_eta2.at(idx_hit) = 1;      // Mark hit as used as rising hit for η2

                        if (dt_time2_trigger < dt_max &&
                            dt_time2_trigger > dt_min) {
                            eta2_hits_counter[layer]++;
                            if (eta2_hit_flags[layer] == false) {
                                eta2_hit_flags[layer] = true;
                            }
                        }
                    }

                    // --------------------------------------------------------------------------------------------
                    // Rising dt calculation subsection (only when both rising time1 and time2 are present)
                    if (time1 != -1 && time2 != -1) {
                        int dt_time1_time2 = time1 - time2;
                        proc_dt_time1_time2.push_back(dt_time1_time2);
                    }

                    // --------------------------------------------------------------------------------------------
                    // Time over Threshold (ToT) calculation subsecion
                    // Step 0: Define time of rising hit in units of ticks
                    int time_rising1 = time1;
                    int time_rising2 = time2;
                    int idx_rising = idx_hit;   // Purely for clarity, since idx_hit already corresponds to the index of the rising hit

                    // Step 1: Look for corresponding falling hits on both sides (same channel, hit_time* != 0, hit is not in hit_used*) and extract time of falling hit
                    int time_falling1 = -1;
                    int time_falling2 = -1;
                    for (size_t idx_falling = idx_rising + 1; idx_falling < n_hits; idx_falling++) {
                        if (hit_channel.at(idx_falling) == hit_channel.at(idx_rising) && hit_rise.at(idx_falling) == 0) {
                            if (hit_time1.at(idx_falling) != 0 &&
                                hit_used_tot_eta1.at(idx_falling) == -1) {

                                time_falling1 = hit_bcid.at(idx_falling) * 30 + hit_time1.at(idx_falling) - 1;
                                hit_used_tot_eta1.at(idx_falling) = 0; // Mark hit as used as falling hit for η1
                            }
                            if (hit_time2.at(idx_falling) != 0 &&
                                hit_used_tot_eta2.at(idx_falling) == -1) {

                                time_falling2 = hit_bcid.at(idx_falling) * 30 + hit_time2.at(idx_falling) - 1;
                                hit_used_tot_eta2.at(idx_falling) = 0; // Mark hit as used as falling hit for η2
                            }
                            break; // Assuming that the first falling hit after the rising hit corresponds to the same signal (which should be the case if the readout window is not too crowded)
                        }
                    } // End of the loop looking for falling hits

                    // Step 2: If corresponding falling hit is found, calculate ToT as time_falling - time_rising
                    if (time_falling1 != -1 && time_rising1 != -1) {
                        int tot1 = time_falling1 - time_rising1;
                        if (tot1 < 0) { // Assuming a positive ToT
                            tot1 = -1; // Set ToT back to -1 if negative (invalid value)
                        }
                        proc_tot1.push_back(tot1);
                    }
                    if (time_falling2 != -1 && time_rising2 != -1) {
                        int tot2 = time_falling2 - time_rising2;
                        if (tot2 < 0) { // Assuming a positive ToT
                            tot2 = -1; // Set ToT back to -1 if negative (invalid value)
                        }
                        proc_tot2.push_back(tot2);
                    }

                    // --------------------------------------------------------------------------------------------
                    // Clusterization subsection

                    // TODO: Move clusterization logic to the front to later use the cluster centers in subsequent calculations
                    // TODO: Add IsClusterCenter vector (-1 = not set, 0 = not center, 1 = cluster center) to later use in the track reconstruction
                    // TODO: In the clusterization logic only do clusterization and calculate the ToT in it's own section

                    // For side η1
                    std::cout << "------------------------------------------------------------------" << std::endl;
                    std::cout << "Processing: Hit: " << idx_hit << "; Layer " << layer << "; Strip: " << strip << "; Time " << time1 << std::endl;

                    // Step 0: Initialize a new cluster containing only the current hit if the hit is not already in a cluster and time information is present
                    if (v_cluster_hit_assinged_eta1.at(idx_hit) == 0 && time1 != -1) {
                        v_cluster_strip_eta1.clear();                       // Clear cluster_strip* vector for the new cluster
                        v_cluster_time_eta1.clear();                        // Clear cluster_time* vector for the new cluster
                        v_cluster_strip_eta1.push_back(strip);              // Push back the strip number of the current hit for later checks of strip adjacency between cluster members
                        v_cluster_time_eta1.push_back(time1);               // Push back the time of the current hit for later time window checks between cluster members and potential cluster partner hits
                        v_cluster_hit_assinged_eta1.at(idx_hit) = 1;              // Mark current hit as belonging to a cluster
                        cluster_hit_used_tot_eta1.at(idx_hit) = 1;            // Mark hit as used as rising hit for η1 in clusterization (to avoid using the same hit again as a rising hit when looking for cluster partners for other hits in ToT calculation)
                        int idx_hit_cluster_center = idx_hit;         // Define current hit as initial cluster center (to be updated later if another hit in the cluster has smaller rising time)

                        // Step 1: Iterate through the following hits of the same event and look for the potential cluster partners
                        for (size_t idx_next_hit = idx_hit + 1; idx_next_hit < n_hits; idx_next_hit++) {
                            if (v_cluster_hit_assinged_eta1.at(idx_next_hit) == 0 &&
                                hit_time1.at(idx_next_hit) != -1 &&
                                hit_rise.at(idx_next_hit) == 1 && 
                                hit_channel.at(idx_next_hit) != trig_channel) { // Only consider rising hits on non-trigger channels as potential cluster partners

                                int next_layer = (hit_channel.at(idx_next_hit) % 24) / 8;
                                if (next_layer != layer) {
                                    // std::cout << "Skipping Hit No. " << idx_next_hit << " since it's on another layer: " << layer << " vs. " << next_layer << std::endl;
                                    continue;
                                }

                                // Only consider hits on the same or neighboring strips as potential cluster partners
                                bool strip_adjacent = false;
                                int next_column = hit_channel.at(idx_next_hit) / 24;
                                int next_strip = 8 * next_column + hit_channel.at(idx_next_hit) % 8;
                                for (size_t idx_cluster_hit = 0; idx_cluster_hit < v_cluster_strip_eta1.size(); idx_cluster_hit++) {
                                    if (abs(next_strip - v_cluster_strip_eta1.at(idx_cluster_hit)) <= 1) {
                                        // If at least one hit in the cluster is on the same or neighboring strip as the potential cluster partner hit, we can consider this potential cluster partner hit for further checks
                                        strip_adjacent = true;
                                        break;
                                    }
                                }
                                if (!strip_adjacent) continue; // If none of the hits in the cluster are on the same or neighboring strip as the potential cluster partner hit, skip this potential cluster partner hits
                                std::cout << "  Potential cluster partner: Hit " << idx_next_hit << "; Layer: " << next_layer << "; Strip: " << next_strip << "; Time " << hit_time1.at(idx_next_hit) << std::endl;

                                // Check if the potential cluster partner hit is within a certain time window (e.g. 15 ns = 18 ticks) from any of the hits already belonging to the cluster
                                int next_time1 = hit_time1.at(idx_next_hit) != 0 ? hit_bcid.at(idx_next_hit) * 30 + hit_time1.at(idx_next_hit) - 1 : -1;
                                if (next_time1 != -1) {
                                    for (size_t idx_cluster_hit = 0; idx_cluster_hit < v_cluster_time_eta1.size(); idx_cluster_hit++) {
                                        if (abs(next_time1 - v_cluster_time_eta1.at(idx_cluster_hit)) < 18) {
                                            v_cluster_hit_assinged_eta1.at(idx_next_hit) = 1; // Mark hit as belonging to a cluster
                                            v_cluster_strip_eta1.push_back(next_strip); // Push back the strip number of the new cluster member hit
                                            v_cluster_time_eta1.push_back(next_time1);  // Push back the time of the new cluster member hit
                                            if (next_time1 < time1) idx_hit_cluster_center = idx_next_hit;
                                            break; // No need to check other hits in the v_cluster_time_eta1 vector, since we just need at least one hit in the cluster to be within the time window from the potential cluster partner hit
                                        }
                                    }
                                }

                            }
                        } // End of loop looking for potential cluster partners

                        // Step 2: If no more potential cluster partners are found for the current hit, store cluster size
                        v_cluster_size_eta1.push_back(v_cluster_time_eta1.size());

                        // DEBUG: Print out parameters of the hits belonging to the cluster for the current hit
                        std::cout << "******************************************************************" << std::endl;
                        std::cout << "Resulting cluster: " << std::endl;
                        for (size_t idx_cluster_hit = 0; idx_cluster_hit < v_cluster_time_eta1.size(); idx_cluster_hit++) {
                            std::cout << idx_cluster_hit << ": strip = " << v_cluster_strip_eta1.at(idx_cluster_hit) << ", time = " << v_cluster_time_eta1.at(idx_cluster_hit) << std::endl;
                        }
                        std::cout << "\n";

                        // Step 3: Calculate the ToT for the cluster center
                        int time_rising_cluster_center = hit_time1.at(idx_hit_cluster_center) != 0 ? hit_bcid.at(idx_hit_cluster_center) * 30 + hit_time1.at(idx_hit_cluster_center) - 1 : -1;
                        int time_falling_cluster_center = -1;

                        for (size_t idx_falling = idx_hit_cluster_center + 1; idx_falling < n_hits; idx_falling++) {
                            if (hit_channel.at(idx_falling) == hit_channel.at(idx_hit_cluster_center) && hit_rise.at(idx_falling) == 0) {
                                if (hit_time1.at(idx_falling) != 0 &&
                                    cluster_hit_used_tot_eta1.at(idx_falling) == -1) {

                                    time_falling_cluster_center = hit_bcid.at(idx_falling) * 30 + hit_time1.at(idx_falling) - 1;
                                    cluster_hit_used_tot_eta1.at(idx_falling) = 0; // Mark hit as used as falling hit for η1
                                }
                                break; // Assuming that the first falling hit after the rising hit corresponds to the same signal (which should be the case if the readout window is not too crowded)
                            }
                        }

                        if (time_rising_cluster_center != -1 && time_falling_cluster_center != -1) {
                            int tot_cluster_center = time_falling_cluster_center - time_rising_cluster_center;
                            if (tot_cluster_center < 0) { // Assuming a positive ToT
                                tot_cluster_center = -1; // Set ToT back to -1 if negative (invalid value)
                            }
                            v_cluster_tot1.push_back(tot_cluster_center);
                        }
                    }   // End of finding a cluster and its ToT for the current hit for the side η1


                    // For side η2
                    // Step 0: Initialize a new cluster containing only the current hit if the hit is not already in a cluster and time information is present
                    if (v_cluster_hit_assinged_eta2.at(idx_hit) == 0 && time2 != -1) {
                        v_cluster_strip_eta2.clear();                       // Clear cluster_strip* vector for the new cluster
                        v_cluster_time_eta2.clear();                        // Clear cluster_time* vector for the new cluster
                        v_cluster_strip_eta2.push_back(strip);              // Push back the strip number of the current hit
                        v_cluster_time_eta2.push_back(time2);               // Push back the time of the current hit
                        v_cluster_hit_assinged_eta2.at(idx_hit) = 1;              // Mark current hit as belonging to a cluster
                        cluster_hit_used_tot_eta2.at(idx_hit) = 1;            // Mark current hit as used as rising hit for η2 in clusterization (to avoid using the same hit again as a rising hit when looking for cluster partners for other hits in ToT calculation)
                        int idx_hit_cluster_center = idx_hit;         // Define current hit as initial cluster center (to be updated later if another hit in the cluster has smaller rising time)

                        // Step 1: Iterate through the following hits of the same event and look for the potential cluster partners
                        for (size_t idx_next_hit = idx_hit + 1; idx_next_hit < n_hits; idx_next_hit++) {
                            if (hit_rise.at(idx_next_hit) == 1 && hit_channel.at(idx_next_hit) != trig_channel) { // Only consider rising hits on non-trigger channels as potential cluster partners

                                int next_layer = (hit_channel.at(idx_next_hit) % 24) / 8;
                                if (next_layer != layer) continue; // Only consider hits in the same layer as potential cluster partners

                                // Only consider hits on the same or neighboring strips as potential cluster partners
                                bool strip_adjacent = false;
                                int next_column = hit_channel.at(idx_next_hit) / 24;
                                int next_strip = 8 * next_column + hit_channel.at(idx_next_hit) % 8;
                                for (size_t idx_cluster_hit = 0; idx_cluster_hit < v_cluster_strip_eta2.size(); idx_cluster_hit++) {
                                    if (abs(next_strip - v_cluster_strip_eta2.at(idx_cluster_hit)) <= 1) {
                                        // If at least one hit in the cluster is on the same or neighboring strip as the potential cluster partner hit, we can consider this potential cluster partner hit for further checks
                                        strip_adjacent = true;
                                        break;
                                    }
                                }
                                if (!strip_adjacent) continue; // If none of the hits in the cluster are on the same or neighboring strip as the potential cluster partner hit, skip this potential cluster partner hits

                                // Check if the potential cluster partner hit is within a certain time window (e.g. 15 ns = 18 ticks) from any of the hits already belonging to the cluster
                                int next_time2 = hit_time2.at(idx_next_hit) != 0 ? hit_bcid.at(idx_next_hit) * 30 + hit_time2.at(idx_next_hit) - 1 : -1;
                                if (v_cluster_hit_assinged_eta2.at(idx_next_hit) == 0 && next_time2 != -1) {
                                    for (size_t idx_cluster_hit = 0; idx_cluster_hit < v_cluster_time_eta2.size(); idx_cluster_hit++) {
                                        if (abs(next_time2 - v_cluster_time_eta2.at(idx_cluster_hit)) < 18) {
                                            v_cluster_hit_assinged_eta2.at(idx_next_hit) = 1; // Mark hit as belonging to a cluster
                                            v_cluster_strip_eta2.push_back(next_strip); // Push back the strip number of the new cluster member hit to cluster_strip* vector for later strip window checks with other potential cluster partners
                                            v_cluster_time_eta2.push_back(next_time2);  // Push back the time of the new cluster member hit to cluster_time* vector for later time window checks with other potential cluster partners
                                            if (next_time2 < time2) idx_hit_cluster_center = idx_next_hit;
                                            break; // No need to check other hits in the v_cluster_time_eta2 vector, since we just need at least one hit in the cluster to be within the time window from the potential cluster partner hit
                                        }
                                    }
                                }

                            }
                        } // End of loop looking for potential cluster partners

                        // Step 2: If no more potential cluster partners are found for the current hit, store cluster size
                        v_cluster_size_eta2.push_back(v_cluster_time_eta2.size());

                        // Step 3: Calculate the ToT for the cluster center
                        int time_rising_cluster_center = hit_time2.at(idx_hit_cluster_center) != 0 ? hit_bcid.at(idx_hit_cluster_center) * 30 + hit_time2.at(idx_hit_cluster_center) - 1 : -1;
                        int time_falling_cluster_center = -1;

                        for (size_t idx_falling = idx_hit_cluster_center + 1; idx_falling < n_hits; idx_falling++) {
                            if (hit_channel.at(idx_falling) == hit_channel.at(idx_hit_cluster_center) && hit_rise.at(idx_falling) == 0) {
                                if (hit_time2.at(idx_falling) != 0 &&
                                    cluster_hit_used_tot_eta2.at(idx_falling) == -1) {

                                    time_falling_cluster_center = hit_bcid.at(idx_falling) * 30 + hit_time2.at(idx_falling) - 1;
                                    hit_used_tot_eta2.at(idx_falling) = 0; // Mark hit as used as falling hit for η2
                                }
                                break; // Assuming that the first falling hit after the rising hit corresponds to the same signal (which should be the case if the readout window is not too crowded)
                            }
                        }

                        if (time_rising_cluster_center != -1 && time_falling_cluster_center != -1) {
                            int tot_cluster_center = time_falling_cluster_center - time_rising_cluster_center;
                            if (tot_cluster_center < 0) { // Assuming a positive ToT
                                tot_cluster_center = -1; // Set ToT back to -1 if negative (invalid value)
                            }
                            v_cluster_tot2.push_back(tot_cluster_center);
                        }
                    }   // End of finding a cluster and its ToT for the current hit for the side η2

                    // --------------------------------------------------------------------------------------------
                    // WIP: Track reconstruction subsection

                } // End of if statement checking for rising hits on non-trigger channel

            } // End of hits loop

            proc_trigger_time.push_back(trigger_time);

            // ----------------------------------------------------------------------------------------------------
            // Efficiency calculation based on hit flags for each layer and update of efficiency counters
            if (trigger_time != -1) { // Only consider events with a valid trigger time (i.e. with a trigger hit)

                // External trigger only efficiency counting
                triggered_events_external++;
                for (int layer = 0; layer < 3; layer++) {
                    if (eta1_hit_flags[layer]) eta1_efficiency_counter[layer]++;
                    if (eta2_hit_flags[layer]) eta2_efficiency_counter[layer]++;
                    if (eta1_hit_flags[layer] || eta2_hit_flags[layer]) eta_or_efficiency_counter[layer]++;
                    if (eta1_hit_flags[layer] && eta2_hit_flags[layer]) eta_and_efficiency_counter[layer]++;
                }

                // External trigger + RPC as a trigger efficiency counting
                // L0
                if ((eta1_hit_flags[1] && eta1_hit_flags[2]) || (eta2_hit_flags[1] && eta2_hit_flags[2])) {
                    triggered_events_rpc[0]++;
                    if (eta1_hit_flags[0]) eta1_efficiency_counter_rpc[0]++;
                    if (eta2_hit_flags[0]) eta2_efficiency_counter_rpc[0]++;
                    if (eta1_hit_flags[0] || eta2_hit_flags[0]) eta_or_efficiency_counter_rpc[0]++;
                    if (eta1_hit_flags[0] && eta2_hit_flags[0]) eta_and_efficiency_counter_rpc[0]++;
                }

                // L1
                if ((eta1_hit_flags[0] && eta1_hit_flags[2]) || (eta2_hit_flags[0] && eta2_hit_flags[2])) {
                    triggered_events_rpc[1]++;
                    if (eta1_hit_flags[1]) eta1_efficiency_counter_rpc[1]++;
                    if (eta2_hit_flags[1]) eta2_efficiency_counter_rpc[1]++;
                    if (eta1_hit_flags[1] || eta2_hit_flags[1]) eta_or_efficiency_counter_rpc[1]++;
                    if (eta1_hit_flags[1] && eta2_hit_flags[1]) eta_and_efficiency_counter_rpc[1]++;
                }

                // L2
                if ((eta1_hit_flags[0] && eta1_hit_flags[1]) || (eta2_hit_flags[0] && eta2_hit_flags[1])) {
                    triggered_events_rpc[2]++;
                    if (eta1_hit_flags[2]) eta1_efficiency_counter_rpc[2]++;
                    if (eta2_hit_flags[2]) eta2_efficiency_counter_rpc[2]++;
                    if (eta1_hit_flags[2] || eta2_hit_flags[2]) eta_or_efficiency_counter_rpc[2]++;
                    if (eta1_hit_flags[2] && eta2_hit_flags[2]) eta_and_efficiency_counter_rpc[2]++;
                }
            } // End of if statement checking for valid trigger time in efficiency counting

            // ----------------------------------------------------------------------------------------------------
            // Fill trees with hit and event information
            input_data_tree->Fill();
            processed_data_tree->Fill();
            clusterization_tree->Fill();

            // Clear hit vectors for next event
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

            v_cluster_size_eta1.clear();
            v_cluster_size_eta2.clear();
            v_cluster_tot1.clear();
            v_cluster_tot2.clear();

            // Reset variables for next event
            n_hits = 0;
            n_clusters_eta1 = 0;
            n_clusters_eta2 = 0;
            n_tracks_eta1 = 0;
            n_tracks_eta2 = 0;
            trigger_time = -1;

        } // End of event processing sections (when clk == 127)

    } // End of main loop processing raw word data


    // ============================================================================================================
    // Efficiency calculation section

    // External trigger only efficiency results
    double eta1_efficiency_external[3] = {0.0, 0.0, 0.0};
    double eta2_efficiency_external[3] = {0.0, 0.0, 0.0};
    double eta_or_efficiency_external[3] = {0.0, 0.0, 0.0};
    double eta_and_efficiency_external[3] = {0.0, 0.0, 0.0};
    
    double eta1_efficiency_external_error[3] = {0.0, 0.0, 0.0};
    double eta2_efficiency_external_error[3] = {0.0, 0.0, 0.0};
    double eta_or_efficiency_external_error[3] = {0.0, 0.0, 0.0};
    double eta_and_efficiency_external_error[3] = {0.0, 0.0, 0.0};

    if (triggered_events_external > 0) {
        for (int layer = 0; layer < 3; layer++) {
            eta1_efficiency_external[layer] = static_cast<double>(eta1_efficiency_counter[layer]) / triggered_events_external;
            eta2_efficiency_external[layer] = static_cast<double>(eta2_efficiency_counter[layer]) / triggered_events_external;
            eta_or_efficiency_external[layer] = static_cast<double>(eta_or_efficiency_counter[layer]) / triggered_events_external;
            eta_and_efficiency_external[layer] = static_cast<double>(eta_and_efficiency_counter[layer]) / triggered_events_external;

            eta1_efficiency_external_error[layer] = sqrt(eta1_efficiency_external[layer] * (1 - eta1_efficiency_external[layer]) / triggered_events_external);
            eta2_efficiency_external_error[layer] = sqrt(eta2_efficiency_external[layer] * (1 - eta2_efficiency_external[layer]) / triggered_events_external);
            eta_or_efficiency_external_error[layer] = sqrt(eta_or_efficiency_external[layer] * (1 - eta_or_efficiency_external[layer]) / triggered_events_external);
            eta_and_efficiency_external_error[layer] = sqrt(eta_and_efficiency_external[layer] * (1 - eta_and_efficiency_external[layer]) / triggered_events_external);
        }
    }

    // External trigger + RPC as a trigger efficiency results
    double eta1_efficiency_rpc[3] = {0.0, 0.0, 0.0};
    double eta2_efficiency_rpc[3] = {0.0, 0.0, 0.0};
    double eta_or_efficiency_rpc[3] = {0.0, 0.0, 0.0};
    double eta_and_efficiency_rpc[3] = {0.0, 0.0, 0.0};

    double eta1_efficiency_rpc_error[3] = {0.0, 0.0, 0.0};
    double eta2_efficiency_rpc_error[3] = {0.0, 0.0, 0.0};
    double eta_or_efficiency_rpc_error[3] = {0.0, 0.0, 0.0};
    double eta_and_efficiency_rpc_error[3] = {0.0, 0.0, 0.0};

    for (int layer = 0; layer < 3; layer++) {
        if (triggered_events_rpc[layer] > 0) {
            eta1_efficiency_rpc[layer] = static_cast<double>(eta1_efficiency_counter_rpc[layer]) / triggered_events_rpc[layer];
            eta2_efficiency_rpc[layer] = static_cast<double>(eta2_efficiency_counter_rpc[layer]) / triggered_events_rpc[layer];
            eta_or_efficiency_rpc[layer] = static_cast<double>(eta_or_efficiency_counter_rpc[layer]) / triggered_events_rpc[layer];
            eta_and_efficiency_rpc[layer] = static_cast<double>(eta_and_efficiency_counter_rpc[layer]) / triggered_events_rpc[layer];

            eta1_efficiency_rpc_error[layer] = sqrt(eta1_efficiency_rpc[layer] * (1 - eta1_efficiency_rpc[layer]) / triggered_events_rpc[layer]);
            eta2_efficiency_rpc_error[layer] = sqrt(eta2_efficiency_rpc[layer] * (1 - eta2_efficiency_rpc[layer]) / triggered_events_rpc[layer]);
            eta_or_efficiency_rpc_error[layer] = sqrt(eta_or_efficiency_rpc[layer] * (1 - eta_or_efficiency_rpc[layer]) / triggered_events_rpc[layer]);
            eta_and_efficiency_rpc_error[layer] = sqrt(eta_and_efficiency_rpc[layer] * (1 - eta_and_efficiency_rpc[layer]) / triggered_events_rpc[layer]);
        }
    }

    // Create histograms to store efficiency results in efficiencies_histograms directory in output ROOT file
    output_file->mkdir("efficiencies_histograms/external_trigger");
    output_file->cd("efficiencies_histograms/external_trigger");

    // External trigger only first
    for (int layer = 0; layer < 3; layer++) {
        TH1F* hist_eff_external_trigger = new TH1F(Form("eff_external_trigger%d", layer), Form("Layer %d Efficiency (External Trigger)", layer), 4, 0.5, 4.5);
        hist_eff_external_trigger->GetXaxis()->SetBinLabel(1, "eta1");
        hist_eff_external_trigger->GetXaxis()->SetBinLabel(2, "eta2");
        hist_eff_external_trigger->GetXaxis()->SetBinLabel(3, "OR");
        hist_eff_external_trigger->GetXaxis()->SetBinLabel(4, "AND");

        hist_eff_external_trigger->SetBinContent(1, eta1_efficiency_external[layer]);
        hist_eff_external_trigger->SetBinError(1, eta1_efficiency_external_error[layer]);
        hist_eff_external_trigger->SetBinContent(2, eta2_efficiency_external[layer]);
        hist_eff_external_trigger->SetBinError(2, eta2_efficiency_external_error[layer]);
        hist_eff_external_trigger->SetBinContent(3, eta_or_efficiency_external[layer]);
        hist_eff_external_trigger->SetBinError(3, eta_or_efficiency_external_error[layer]);
        hist_eff_external_trigger->SetBinContent(4, eta_and_efficiency_external[layer]);
        hist_eff_external_trigger->SetBinError(4, eta_and_efficiency_external_error[layer]);
    }

    // Then external trigger + RPC as a trigger
    output_file->cd();
    output_file->mkdir("efficiencies_histograms/external_plus_rpc_trigger");
    output_file->cd("efficiencies_histograms/external_plus_rpc_trigger");

    for (int layer = 0; layer < 3; layer++) {
        TH1F* hist_eff_rpc_trigger = new TH1F(Form("eff_rpc_layer%d", layer), Form("Layer %d Efficiency (External and RPC Trigger)", layer), 4, 0.5, 4.5);
        hist_eff_rpc_trigger->GetXaxis()->SetBinLabel(1, "eta1");
        hist_eff_rpc_trigger->GetXaxis()->SetBinLabel(2, "eta2");
        hist_eff_rpc_trigger->GetXaxis()->SetBinLabel(3, "OR");
        hist_eff_rpc_trigger->GetXaxis()->SetBinLabel(4, "AND");

        hist_eff_rpc_trigger->SetBinContent(1, eta1_efficiency_rpc[layer]);
        hist_eff_rpc_trigger->SetBinError(1, eta1_efficiency_rpc_error[layer]);
        hist_eff_rpc_trigger->SetBinContent(2, eta2_efficiency_rpc[layer]);
        hist_eff_rpc_trigger->SetBinError(2, eta2_efficiency_rpc_error[layer]);
        hist_eff_rpc_trigger->SetBinContent(3, eta_or_efficiency_rpc[layer]);
        hist_eff_rpc_trigger->SetBinError(3, eta_or_efficiency_rpc_error[layer]);
        hist_eff_rpc_trigger->SetBinContent(4, eta_and_efficiency_rpc[layer]);
        hist_eff_rpc_trigger->SetBinError(4, eta_and_efficiency_rpc_error[layer]);
    }

    // ============================================================================================================
    // Write output and cleanup
    output_file->Write();
    output_file->Close();

    tmp_file.close();
}