#include "dataAnalyzer.hpp"


int main(int argc, char** argv) {

    // Process command line arguments
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <dt_max> <dt_min> <--self>" << std::endl;
        return 1;
    }
    std::string input = argv[1];
    int dt_max, dt_min;

    // Throw an error if input is a directory
    if (std::filesystem::exists(input) && std::filesystem::is_directory(input)) {
        std::cerr << "Error: Input path is a directory. Please provide a single file as input." << std::endl;
        return 1;
    }
    
    try {
        dt_max = std::stoi(argv[2]);
        dt_min = std::stoi(argv[3]);
    } catch (const std::exception& e) {
        std::cerr << "Error: dt_max and dt_min must be valid integers" << std::endl;
        return 1;
    }
    
    bool self_trigger = (argc > 4 && std::string(argv[4]) == "--self");
    
    // Suppress unused variable warnings (TODO: integrate these parameters into analysis)
    (void)self_trigger;

    // Setup output file and trees
    DataAnalyzer analyzer;
    analyzer.setupOutputFile();
    analyzer.setupBranches();

    // Process input file
    try {
        analyzer.processInputData(input, dt_max, dt_min);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error processing input data: " << e.what() << std::endl;
        return 1;
    }
}