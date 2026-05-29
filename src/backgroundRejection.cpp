#include "backgroundRejection.hpp"


TH1* readDataAsHist(TTree* tree, TBranch* branch) {
    // Check if tree and branch are valid
    if (!tree || !branch) {
        std::cerr << "Error: Invalid tree or branch pointer." << std::endl;
        return nullptr;
    }

    // Set branch address
    int value;
    branch->SetAddress(&value);

    // Find min/max
    int min = std::numeric_limits<int>::max();
    int max = std::numeric_limits<int>::lowest();
    Long64_t nEntries = branch->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (value < min) min = value;
        if (value > max) max = value;
    }

    // Create histogram and fill
    TH1* hist = new TH1D("hist", "Data Histogram", max - min + 1, min - 0.5, max + 0.5);
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        hist->Fill(value);
    }

    return hist;
}


TF1* flatPlusGaussian() {
    // Create TF1 with formula string
    // Parameters: [0]=A_flat, [1]=A_sig, [2]=mu, [3]=sigma
    TF1* f1 = new TF1("flatPlusGaussian", "[0] + [1]*exp(-0.5*((x-[2])/[3])^2)");
    
    // Set parameter names
    f1->SetParNames("A_flat", "A_sig", "mu", "sigma");
    
    return f1;
}


std::vector<double> fitBackground(TH1* hist, TF1* fitFunc, const std::vector<double> initialParams) {
    // Check if histogram and fit function are valid
    if (!hist || !fitFunc) {
        std::cerr << "Error: Invalid histogram or fit function pointer." << std::endl;
        return {};
    }

    // Set initial parameters for the fit function
    for (size_t i = 0; i < initialParams.size(); ++i) {
        fitFunc->SetParameter(i, initialParams[i]);
    }

    // Fit the histogram with the provided function
    TFitResultPtr result = hist->Fit(fitFunc, "QS");
    if (result->Status() != 0) {
        std::cerr << "Warning: Fit may not have converged properly. Status: " << result->Status() << std::endl;
    }

    // Extract fitted parameters
    std::vector<double> fittedParams;
    for (size_t i = 0; i < initialParams.size(); ++i) {
        fittedParams.push_back(fitFunc->GetParameter(i));
    }

    return fittedParams;
}


std::pair<int, int> setSignalRange(TF1* fitFunc) {
    // Check if fit function is valid    
    if (!fitFunc) {
        std::cerr << "Error: Invalid fit function pointer." << std::endl;
        return {0, 0};
    }

    double mu = fitFunc->GetParameter("mu");
    double sigma = fitFunc->GetParameter("sigma");

    int lowerBound = static_cast<int>(mu - 3 * sigma);
    int upperBound = static_cast<int>(mu + 3 * sigma);
    return std::make_pair(lowerBound, upperBound);
}


struct SelectionMask {
    std::string name;
    std::vector<bool> mask;

    SelectionMask(const std::string& n, std::vector<bool>&& m) 
        : name(n), mask(std::move(m)) {}
};


std::vector<SelectionMask> getMasks(
    const std::vector<int>& hit_time1,
    const std::vector<int>& hit_time2,
    const std::vector<int>& hit_rise) {

    std::vector<SelectionMask> masks;
    size_t nHits = hit_time1.size();

    // Add eta2_rising
    masks.push_back({"eta2_rising", std::vector<bool>(nHits)});
    for (size_t i = 0; i < nHits; ++i) {
        masks.back().mask[i] = (hit_time2[i] > 0) && (hit_rise[i] == 1);
    }

    // Add eta2_falling
    masks.push_back({"eta2_falling", std::vector<bool>(nHits)});
    for (size_t i = 0; i < nHits; ++i) {
        masks.back().mask[i] = (hit_time2[i] > 0) && (hit_rise[i] == 0);
    }

    return masks;
}


void createSignalBranch(TTree* tree, const std::map<std::string, std::pair<double, double>>& categoryRanges) {

    // Set up branches
    std::vector<int>* hit_clk = nullptr;
    std::vector<int>* hit_time1 = nullptr;
    std::vector<int>* hit_time2 = nullptr;
    std::vector<int>* hit_rise = nullptr;

    if (!tree->GetBranch("hit_clk") || !tree->GetBranch("hit_time1") ||
        !tree->GetBranch("hit_time2") || !tree->GetBranch("hit_rise")) {
        std::cerr << "Error: Required branches not found in tree." << std::endl;
        return;
    }

    tree->SetBranchAddress("hit_clk", &hit_clk);
    tree->SetBranchAddress("hit_time1", &hit_time1);
    tree->SetBranchAddress("hit_time2", &hit_time2);
    tree->SetBranchAddress("hit_rise", &hit_rise);

    // Create output branch
    std::vector<int> signalFlags;
    TBranch* outputBranch = tree->Branch("is_signal", &signalFlags);

    Long64_t nEntries = tree->GetEntries();

    for (Long64_t event_idx = 0; event_idx < nEntries; ++event_idx) {
        tree->GetEntry(event_idx);

        auto masks = getMasks(*hit_time1, *hit_time2, *hit_rise);
        size_t nHits = hit_clk->size();
        signalFlags.clear();
        signalFlags.resize(nHits, 0);

        // Apply each selection's range
        for (const auto& sel : masks) {
            auto it = categoryRanges.find(sel.name);
            if (it == categoryRanges.end()) continue; // No range defined for this selection 

            double lo = it->second.first;
            double hi = it->second.second;
            
            for (size_t hit_idx = 0; hit_idx < nHits; ++hit_idx) {
                if (sel.mask[hit_idx] && 
                    (*hit_clk)[hit_idx] >= lo && (*hit_clk)[hit_idx] <= hi) {
                    signalFlags[hit_idx] = 1;
                }
            }
        }

        outputBranch->Fill();
    }
}