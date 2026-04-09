#include "dataAnalyzer.hpp"


// ==========================================================================================
// DataAnalyzer class implementation
// ==========================================================================================
/// Constructor and destructor for DataAnalyzer class
DataAnalyzer::DataAnalyzer() {
    output_file = nullptr;
    input_data_tree = nullptr;
    processed_data_tree = nullptr;
    clusterization_tree = nullptr;
    track_reconstruction_tree = nullptr;
    current_event = nullptr;
    current_event_number = 0;
    n_hits = 0;
    BC0 = 0;
}

DataAnalyzer::~DataAnalyzer() {
    if (output_file) { 
        // Write all trees to file before closing
        if (input_data_tree) input_data_tree->Write();
        if (processed_data_tree) processed_data_tree->Write();
        if (clusterization_tree) clusterization_tree->Write();
        if (track_reconstruction_tree) track_reconstruction_tree->Write();
        output_file->Close(); 
        delete output_file; 
    }
}

/// Utility function to setup the output ROOT file and create the necessary trees
void DataAnalyzer::setupOutputFile() {
    output_file = new TFile("output.root", "RECREATE");
    input_data_tree = new TTree("InputData", "Raw hit data from the DCT");
    processed_data_tree = new TTree("ProcessedData", "Processed hit and event-level data");
    clusterization_tree = new TTree("Clusterization", "Cluster-level data");
    track_reconstruction_tree = new TTree("TrackReconstruction", "Track-level data");
}

/// Utility function to setup branches for all trees in the output file
void DataAnalyzer::setupBranches() {
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
    processed_data_tree->Branch("n_events", &current_event_number);
    processed_data_tree->Branch("n_hits", &n_hits);
    processed_data_tree->Branch("proc_layer", &proc_layer);
    processed_data_tree->Branch("proc_strip", &proc_strip);
    processed_data_tree->Branch("proc_time1", &proc_time1);
    processed_data_tree->Branch("proc_time2", &proc_time2);
    processed_data_tree->Branch("proc_dt_time1_time2", &proc_dt_time1_time2);
    processed_data_tree->Branch("proc_trigger_time", &proc_trigger_time);
    processed_data_tree->Branch("proc_dt_time1_trigger", &proc_dt_time1_trigger);
    processed_data_tree->Branch("proc_dt_time2_trigger", &proc_dt_time2_trigger);
    processed_data_tree->Branch("proc_tot1", &proc_tot1);
    processed_data_tree->Branch("proc_tot2", &proc_tot2);

    // Branch definitions for clusterization tree
    clusterization_tree->Branch("cluster_size_eta1", &cluster_size_eta1);
    clusterization_tree->Branch("cluster_size_eta2", &cluster_size_eta2);
    clusterization_tree->Branch("cluster_tot1", &cluster_tot1);
    clusterization_tree->Branch("cluster_tot2", &cluster_tot2);

    // Branch definitions for track reconstruction tree
    track_reconstruction_tree->Branch("track_length_eta1", &track_length_eta1);
    track_reconstruction_tree->Branch("track_length_eta2", &track_length_eta2);
}

/// Utility function to clear all event-level vectors before processing the next event
void DataAnalyzer::clearEventVectors() {
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

    cluster_size_eta1.clear();
    cluster_size_eta2.clear();
    cluster_tot1.clear();
    cluster_tot2.clear();

    track_length_eta1.clear();
    track_length_eta2.clear();
}

/// Utility function to push raw hit data into the corresponding vectors for tree filling
void DataAnalyzer::pushBackHitData(const Hit& hit) {
    hit_clk.push_back(hit.getClk());
    hit_channel.push_back(hit.getChannel());
    hit_raw_bcid.push_back(hit.getRawBCID());
    hit_bcid.push_back(hit.getBCID());
    hit_time1.push_back(hit.getRawTimeEta1());
    hit_time2.push_back(hit.getRawTimeEta2());
    hit_rise.push_back(hit.getRise());
    hit_raw_bcout.push_back(hit.getRawBCOut());
    hit_bcout.push_back(hit.getBCOut());
}

/// Utility function to push processed hit data into the corresponding vectors for tree filling
void DataAnalyzer::pushBackProcessedData(const Event& event) {

    proc_trigger_time.push_back(event.getTriggerTime());
    // Loop over hits in the event and push processed data into vectors
    // Note: Only push hits with rise == 1 (to match with ToT calculation which only processes rising edges)
    for (const auto& hit : event.getHits()) {
        // Only process non-trigger channel and rising edges for consistency with ToT calculation
        if (hit.getChannel() == event.getTriggerChannel() || hit.getRise() != 1) continue;
        
        if (hit.hasEta1Time() || hit.hasEta2Time()) {
            proc_layer.push_back(hit.getLayer());
            proc_strip.push_back(hit.getStrip());
        }
        if (hit.hasEta1Time()) {
            proc_time1.push_back(hit.getTimeEta1());
            proc_dt_time1_trigger.push_back(hit.getTimeEta1() - event.getTriggerTime());
        }
        if (hit.hasEta2Time()) {
            proc_time2.push_back(hit.getTimeEta2());
            proc_dt_time2_trigger.push_back(hit.getTimeEta2() - event.getTriggerTime());
        }
        if (hit.hasEta1Time() && hit.hasEta2Time()) proc_dt_time1_time2.push_back(hit.getTimeEta1() - hit.getTimeEta2());
        if (hit.getToT1() > 0) proc_tot1.push_back(hit.getToT1());
        if (hit.getToT2() > 0) proc_tot2.push_back(hit.getToT2());
    }
}

/// TODO: Utility function to push cluster-level data into the corresponding vectors for tree filling
void DataAnalyzer::pushBackClusterDataEta1(const Cluster& cluster) {
    cluster_size_eta1.push_back(cluster.getSize()); // Total cluster size
    if (cluster.getTot1() > 0) cluster_tot1.push_back(cluster.getTot1());
}

void DataAnalyzer::pushBackClusterDataEta2(const Cluster& cluster) {
    cluster_size_eta2.push_back(cluster.getSize()); // Total cluster size
    if (cluster.getTot2() > 0) cluster_tot2.push_back(cluster.getTot2());
}

/// TODO: Utility function to push track-level data into the corresponding vectors for tree filling
void DataAnalyzer::pushBackTrackData(const Track& track) {
    track_length_eta1.push_back(track.getLayerCount(0)); // Eta1 layer count
    track_length_eta2.push_back(track.getLayerCount(1)); // Eta2 layer count
}

/// Main entry point for processing input data from file or directory
void DataAnalyzer::processInputData(const std::string& input_path) {
    namespace fs = std::filesystem;

    if (!fs::exists(input_path)) {
        std::cerr << "ERROR: Input path does not exist: " << input_path << std::endl;
        return;
    }

    // Check if input_path is a file or directory
    if (fs::is_regular_file(input_path)) {
        // Single file - process it directly
        std::cout << "Processing single file: " << input_path << std::endl;
        processFile(input_path);
    } else if (fs::is_directory(input_path)) {
        // Directory - process all files in it
        std::cout << "Processing files in directory: " << input_path << std::endl;
        for (const auto& entry : fs::directory_iterator(input_path)) {
            if (fs::is_regular_file(entry)) {
                std::cout << "Processing file: " << entry.path() << std::endl;
                processFile(entry.path().string());
            }
        }
    } else {
        std::cerr << "ERROR: Input path is neither a file nor a directory!" << std::endl;
    }

    // Process any remaining event at the end of input
    if (current_event_hits.size() > 0) {
        processEvent();
    }
}

/// Process a single file by reading lines and extracting word data
void DataAnalyzer::processFile(const std::string& file_path) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open file: " << file_path << std::endl;
        return;
    }

    int clk, word, raw_bcout;

    // Read words from file: "clk word raw_bcout" format (all in hexadecimal or as they appear
    while (infile >> std::dec >> clk >> std::hex >> word >> raw_bcout) {
        processSingleWord(clk, word, raw_bcout);
    }

    infile.close();
}

/// Helper function to extract raw BCID from DCT word
int DataAnalyzer::extractRawBCID(int word) {
    int rise = word & 0x01;               // Bit 0
    int raw_bcid;
    
    if (rise == 1) {  // Rising edge
        raw_bcid = (word >> 12) & 0xFF;   // Bits 12-19
    } else {          // Falling edge
        raw_bcid = (word >> 11) & 0x1FF;  // Bits 11-19
    }
    
    return raw_bcid;
}

/// Process a single word and accumulate into events
/// When CLK reaches 127 (end of BC frame), process the accumulated event
void DataAnalyzer::processSingleWord(int clk, int word, int raw_bcout) {
    // Only process not empty words
    if (word != EMPTY_WORD) {
        // Calculate BC0 if this is the first hit of an event
        if (current_event_hits.empty()) {
            int raw_bcid = extractRawBCID(word);
            BC0 = raw_bcid % 256;
        }

        // Create a Hit object from the raw word data
        Hit hit(clk, word, raw_bcout, BC0);
        hit.setIdx(current_event_hits.size());
        current_event_hits.push_back(hit);

        // Accumulate raw hit data into vector
        pushBackHitData(hit);
    }

    // Check for end of event and process it
    if (clk == 127 && current_event_hits.size() > 0) {
        processEvent();
        current_event_number++;
        current_event_hits.clear();

        // Clear processed data vectors for the next event
        clearEventVectors();
    }
}

/// Process a complete event that has been accumulated
void DataAnalyzer::processEvent() {
    n_hits = current_event_hits.size();
    std::cout << "\n" << std::endl;
    std::cout << "╔═════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║ Processing event No. " << current_event_number << " with " << n_hits << " hits" << std::endl;
    
    // Fill the input data tree with raw hit information
    if (input_data_tree && hit_clk.size() > 0) {
        input_data_tree->Fill();
    }

    // Create an Event object and transfer hits using move semantics
    Event event(current_event_number, std::move(current_event_hits));

    // Extract trigger information from the event
    event.extractTriggerTime();

    // Calculate ToT before clusterization
    event.calculateTOT();

    // Push back processed hit data into vectors for tree filling
    pushBackProcessedData(event);

    // Fill the processed data tree
    if (processed_data_tree && proc_layer.size() > 0) {
        processed_data_tree->Fill();
    }

    // Clusterization
    event.clusterize();

    // Calculate ToT for cluster centers
    event.calculateTOTCluster();

    // Push back cluster-level data for both eta1 and eta2 sides
    for (const auto& cluster : event.getClustersEta1()) {
        pushBackClusterDataEta1(cluster);
    }
    for (const auto& cluster : event.getClustersEta2()) {
        pushBackClusterDataEta2(cluster);
    }

    // Fill the clusterization tree (write if either eta1 or eta2 has clusters)
    if (clusterization_tree && (cluster_size_eta1.size() > 0 && cluster_size_eta2.size() > 0)) {
        clusterization_tree->Fill();
    }

    /// WIP: Track reconstruction
    event.reconstructTracks();

    /// TODO: Implement track reconstruction, efficiency calculation
    /// UPDATE: Implement push backs for processed data vectors using the new pushBackProcessedData method
    // - Call current_event->reconstructTracks()
    // - Call current_event->calculateEfficiency()
    // - Call current_event->calculateTOTCluster() for each cluster
    // - Fill processed_data_tree, clusterization_tree, track_reconstruction_tree
}