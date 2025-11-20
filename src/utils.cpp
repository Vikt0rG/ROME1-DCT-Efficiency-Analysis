#include <TFile.h>

#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <vector>

#include "utils.hpp"


std::vector<Int_t> numberOfHitsPerEvent(
    TTree* rawDataTree,
    const TString& hitBranchName,
    bool verbose) 
{
    // Get the total number of entries
    Long64_t nEntries = rawDataTree->GetEntries();
    std::cout << "Number of entries in raw data TTree \"" << rawDataTree->GetName() << "\" : " << nEntries << std::endl;

    // Vector to hold the number of hits per event
    std::vector<Int_t> hitsPerEvent = std::vector<Int_t>(nEntries, -1);
    Int_t nHits = 0;
    rawDataTree->SetBranchAddress(hitBranchName, &nHits);
    for (Int_t i = 0; i < nEntries; ++i) {
        rawDataTree->GetEntry(i);
        hitsPerEvent[i] = nHits;
    }

    // Print out number of hits per event
    if (verbose) {
        for (Long64_t i = 0; i < nEntries; ++i) {
            std::cout << "Event "<< i << ": " << hitsPerEvent[i] << " hits" << std::endl;
        }
    }

    return hitsPerEvent;
}


std::vector<Int_t> mapChannelToLayer(
    TTree* rawDataTree,
    const TString& hitChannelBranchName,
    bool verbose)
{

    // Placeholder output
    return std::vector<Int_t>();
}