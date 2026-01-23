{
    // ROOT macro structure with curly braces

    gROOT->ProcessLine("#include \"src/utils.cpp\"");
    gROOT->ProcessLine("#include \"<vector>\"");

    // TODO: Unused in orignal code - remove if not needed
    // Edit the cuts to be applied to the difference between hit time and trigger time to compute efficiency (values depend on cabling)
    int dt_max = -100;
    int dt_min = -180;

    ifstream tmp_file;
    tmp_file.open("tmp.strip");

    TFile *output_file = TFile::Open("out.root", "RECREATE");
    if (!output_file) return;

    // Global constants
    const int trig_channel = 143;    // Trigger signal from scintillator is treated as a DCT hit and is sent to channel 143

    // Tree
    int event = 0;
    int nHits = 0;

    // Word structure taken from DCT 'documentation' at https://gitlab.cern.ch/groups/atlas-tdaq-muon-barrel-trigger/-/wikis/home/DCT/%7BBI-DCT-test-station-at-BB5%7D
    std::vector<int> hit_channel;    // 8-bit strip ID
    std::vector<int> hit_bcid;       // 8-bit (or 9-bit for falling edge) BC ID
    std::vector<int> hit_time1;      // 5-bit time of η1
    std::vector<int> hit_time2;      // 5-bit time of η2
    std::vector<int> hit_rise;       // 1-bit rise or fall 'status' (1 = rise, 0 = fall)

    std::vector<int> hit_clk;        // Time of hit in readout window (tick = 3.125 ns)
    std::vector<int> hit_bcout;      // Number of BC at which word exits the DCT
}

// TODO: Pack into a function in utils.cpp
// Extraction logic from raw words:
int word = /* raw hex data word from DCT */;
int channel = (word >> 20) & 0xFF;    // Bits 20-27
int rise = word & 0x01;               // Bit 0

if (rise == 1) {                      // Rising edge
    bcid = (word >> 12) & 0xFF;       // Bits 12-19
    time1 = (word >> 7) & 0x1F;       // Bits 7-11
    time2 = (word >> 1) & 0x1F;       // Bits 1-5

} else {                              // Falling edge
    bcid = (word >> 11) & 0x1FF;      // Bits 11-19
    time1 = (word >> 6) & 0x1F;       // Bits 6-10
    time2 = (word >> 1) & 0x1F;       // Bits 1-5
}

