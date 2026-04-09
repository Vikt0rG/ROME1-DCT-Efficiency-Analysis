#include "dataAnalyzer.hpp"


int main(int argc, char** argv) {

    // Get command line arguments
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file/input_directory> <dt_max> <dt_min> <--self>" << std::endl;
        return 1;
    }

    std::string input = argv[1];
    int dt_max, dt_min;
    
    try {
        dt_max = std::stoi(argv[2]);
        dt_min = std::stoi(argv[3]);
    } catch (const std::exception& e) {
        std::cerr << "Error: dt_max and dt_min must be valid integers" << std::endl;
        return 1;
    }
    
    bool self_trigger = (argc > 4 && std::string(argv[4]) == "--self");
    
    // Suppress unused variable warnings (TODO: integrate these parameters into analysis)
    (void)dt_max;
    (void)dt_min;
    (void)self_trigger;

    // Setup output file and trees
    DataAnalyzer analyzer;
    analyzer.setupOutputFile();
    analyzer.setupBranches();

    // Process input data from file or directory
    analyzer.processInputData(input);

}