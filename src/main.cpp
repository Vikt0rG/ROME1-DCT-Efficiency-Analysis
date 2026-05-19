#include "dataProcesser.hpp"


int main(int argc, char** argv) {

    // Process command line arguments
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <dt_max> <dt_min> [--self] [--no-external] [--use-old-data] [--output=<output_file>]" << std::endl;
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

    bool self_trigger = false;
    bool use_external_trigger = true;
    DataProcesser::InputFormat input_format = DataProcesser::InputFormat::FiledumpPackets;

    for (int i = 4; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--self") {
            self_trigger = true;
        } else if (arg == "--no-external") {
            use_external_trigger = false;
        } else if (arg == "--use-old-data") {
            input_format = DataProcesser::InputFormat::DecodedWords;
        }
    }

    // Suppress unused variable warnings (TODO: integrate these parameters into analysis)
    (void)self_trigger;

    // Setup output file and trees
    DataProcesser processor;
    processor.setupOutputFile();
    processor.setupBranches();

    // Process input file
    try {
        processor.processInputData(input, dt_max, dt_min, input_format, use_external_trigger);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error processing input data: " << e.what() << std::endl;
        return 1;
    }

    // Analyze data and produce plots
}