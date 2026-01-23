#include <iostream>
#include <TFile.h>
#include <TTree.h>

#include <utils.hpp>


// Main function is there purely to test utility functions
// Actual idea is to have ROOT macros in root_macros/ for event processing

int main(int argc, char** argv) {
    // Arguments:
    // argv[1] : input ROOT file name
    // argv[2] : option to print out intermediate results [bool]


    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input ROOT file>" " <verbose (0 or 1)>" << std::endl;
        return 1;
    }

    // Parse arguments
    TString inputFileName = argv[1];
    bool verbose = static_cast<bool>(std::stoi(argv[2]));

    // Open the input ROOT file and verify tree and branch names
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
    }

    TString rawDataTreeName = "tree"; // Assuming the TTree name is "tree"
    TTree* rawDataTree = dynamic_cast<TTree*>(inputFile->Get(rawDataTreeName));
    if (!rawDataTree) {
        std::cerr << "Error: TTree '" << rawDataTreeName << "' not found in file " << inputFileName << std::endl;
        inputFile->Close();
    }

    TString hitBranchName = "nHits";  // Assuming the hit branch name is "nHits"
    if (!rawDataTree->GetBranch(hitBranchName)) {
        std::cerr << "Error: Branch '" << hitBranchName << "' not found " << std::endl;
        inputFile->Close();
    }

    if (verbose) std::cout << std::string(40, '-') << std::endl;
    numberOfHitsPerEvent(rawDataTree, hitBranchName, verbose);
    if (verbose) std::cout << std::string(40, '-') << std::endl;

    return 0;
}