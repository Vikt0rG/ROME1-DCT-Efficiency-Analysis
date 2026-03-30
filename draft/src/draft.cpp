#include <draft.hpp>

// ==========================================================================================
// Hit class implementation
// ==========================================================================================
/// Constructor that takes raw DCT word and BC0 for processing
Hit::Hit(int clk, int word, int BC0) : clk(clk) {
    // Decode the raw DCT word to extract channel, BCID, time information and rise/fall edge
    decodeDCTWord(word);

    // Apply BC0 in the event processing
    applyBCWrapAround(BC0);
    
    // Convert raw to physical time
    time_eta1 = TimeUtils::convertRawTimeToPhysical(raw_time_eta1, bcid);
    time_eta2 = TimeUtils::convertRawTimeToPhysical(raw_time_eta2, bcid);

    // Map channel to layer and strip
    geometryMapping();
}

/// Utility function for extracting information from the raw DCT word
void Hit::decodeDCTWord(int word) {
    channel = (word >> 20) & 0xFF;          // Bits 20-27
    rise = word & 0x01;                     // Bit 0

    if (rise == 1) {                      // Rising edge
        raw_bcid = word >> 12 & 0xFF;       // Bits 12-19
        raw_time_eta1 = word >> 7 & 0x1F;   // Bits 7-11
        raw_time_eta2 = word >> 1 & 0x3F;   // Bits 1-5

    } else {                              // Falling edge
        raw_bcid = word >> 11 & 0x1FF;      // Bits 11-19
        raw_time_eta1 = word >> 6 & 0x1F;   // Bits 6-10
        raw_time_eta2 = word >> 1 & 0x1F;   // Bits 1-5
    }
}

/// Utility function for applying BC wrap-around correction to BCID and BCOut
void Hit::applyBCWrapAround(int BC0) {
    bcid = raw_bcid - BC0;
    if (bcid < -128) bcid += 256;
    if (bcid > 128) bcid -= 256;
}

/// Utility function for mapping channel number to detector layer and strip
void Hit::geometryMapping() {
    int layer = (channel % 24) / 8;         // Detector layer (0, 1 or 2)
    int column = channel / 24;              // Detector column (0, 1, 2, ... 5)
    int strip = 8 * column + channel % 8;   // Strip number (0-47)

    this->layer = layer;
    this->strip = strip;
}

// ==========================================================================================
// Event class implementation
// ==========================================================================================
/// Constructor that initializes the event number and sets default values for trigger information
Event::Event(int event_number) : event_number(event_number) {
    trigger_time = -1;      // Initialize trigger time to invalid
    trigger_channel = 143;  // Initialize trigger channel
}

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
}

DataAnalyzer::~DataAnalyzer() {
    if (output_file) { output_file->Close(); delete output_file; }
    if (input_data_tree) { delete input_data_tree; }
    if (processed_data_tree) { delete processed_data_tree; }
    if (clusterization_tree) { delete clusterization_tree; }
    if (track_reconstruction_tree) { delete track_reconstruction_tree; }
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
