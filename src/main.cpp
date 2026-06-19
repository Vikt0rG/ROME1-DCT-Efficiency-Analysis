#include "dataProcesser.hpp"
#include "dataAnalyzer.hpp"
#include "dataPlotter.hpp"


int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " process <input_file> <dt_max> <dt_min> [--no-external] [--use-old-data] [--reject-background] [--self]\n"
                  << "  " << argv[0] << " analyze --config <config_file>" << std::endl
                  << "  " << argv[0] << " plotter --config <config_file>" << std::endl;
        return 1;
    }

    std::string command = argv[1];

    // Command 1: Process input data and produce output ROOT file with trees
    if (command == "process") {
        if (argc < 5) {
            std::cerr << "Usage: " << argv[0] << " process <input_file> <dt_max> <dt_min> [--no-external] [--use-old-data] [--reject-background] [--self]" << std::endl;
            return 1;
        }

        std::string input = argv[2];
        int dt_max = 0;
        int dt_min = 0;

        if (std::filesystem::exists(input) && std::filesystem::is_directory(input)) {
            std::cerr << "Error: Input path is a directory. Please provide a single file as input." << std::endl;
            return 1;
        }

        try {
            dt_max = std::stoi(argv[3]);
            dt_min = std::stoi(argv[4]);
        } catch (const std::exception&) {
            std::cerr << "Error: dt_max and dt_min must be valid integers" << std::endl;
            return 1;
        }

        bool use_external_trigger = true; // a.k.a. no_external = false
        bool reject_background = false; // By default, apply background rejection to improve efficiency results
        bool self_trigger = false;
        DataProcesser::InputFormat input_format = DataProcesser::InputFormat::FiledumpPackets;

        for (int i = 5; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--no-external") {
                use_external_trigger = false;
            } else if (arg == "--use-old-data") {
                input_format = DataProcesser::InputFormat::DecodedWords;
            } else if (arg == "--reject-background") {
                reject_background = true;
            } else if (arg == "--self") {
                self_trigger = true;
            }
        }

        // Suppress unused variable warnings (TODO: integrate these parameters into analysis)
        (void)self_trigger;

        DataProcesser processor(input, dt_max, dt_min, input_format, use_external_trigger, reject_background);
        processor.setupOutputFile();
        processor.setupBranches();

        try {
            processor.processInputData();
        } catch (const std::exception& e) {
            std::cerr << "Error processing input data: " << e.what() << std::endl;
            return 1;
        }

        return 0;
    }

    // Command 2: Analyze processed data and produce summary statistics
    if (command == "analyze") {
        std::string config_file;
        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--config" && i + 1 < argc) {
                config_file = argv[++i];
            }
        }

        if (config_file.empty()) {
            std::cerr << "Usage: " << argv[0] << " analyze --config <config_file>" << std::endl;
            return 1;
        }

        try {
            DataAnalyzer analyzer(config_file);
            analyzer.produceSummaryStats();
        } catch (const std::exception& e) {
            std::cerr << "Error analyzing data: " << e.what() << std::endl;
            return 1;
        }

        return 0;
    }

    // Command 3: Produce summary plots from processed data
    if (command == "plotter") {
        std::vector<std::string> config_paths;
        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--config" && i + 1 < argc) {
                config_paths.push_back(argv[++i]);
            }
        }

        if (config_paths.empty()) {
            std::cerr << "Usage: " << argv[0] << " plotter --config <config_files>" << std::endl;
            return 1;
        }

        try {
            DataPlotter plotter(config_paths);
            plotter.produceSummaryPlots();
        } catch (const std::exception& e) {
            std::cerr << "Error plotting data: " << e.what() << std::endl;
            return 1;
        }

        return 0;
    }

    std::cerr << "Unknown command: " << command << "\n"
              << "Valid commands are: process, analyze, plotter" << std::endl;
    return 1;
}