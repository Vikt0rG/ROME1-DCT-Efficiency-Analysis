#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>

#include <vector>
#include <dirent.h>
#include <regex>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <cmath>

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include "TLatex.h"
#include "TLine.h"


// NOTE: Run macro with command like: 
//       root -l -b -q 'root_macros/analysis.cpp("data/2026-01-28")'
//       to run the analysis on all files in the data directory

void analysis (const char* input_data_dir) {

    // ============================================================================================================
    // Efficiency curve and Working Point estimation

    // Open input directory
    std::string input_dir(input_data_dir);
    if(input_dir.empty()) return;
    if(input_dir.back() != '/') input_dir.push_back('/');

    DIR *dir;
    struct dirent *ent;
    dir = opendir(input_dir.c_str());
    if(!dir) {
        std::cerr << "Cannot open directory: " << input_dir << std::endl;
        return;
    }

    // Get data from the directory name (assuming a date in format like YYYY-MM-DD is present in the name)
    std::string date_str;
    std::regex date_re(".*(\\d{4}-\\d{2}-\\d{2}).*");
    std::smatch date_match;
    if(std::regex_search(input_dir, date_match, date_re)) {
        date_str = date_match[1].str();
        std::cout << "Extracted date from input directory name: " << date_str << std::endl;
    } else {
        std::cout << "No date found in input directory name." << std::endl;
    }

    // Filter files for external and self-trigger datasets
    std::vector<std::string> external_trigger_files;
    std::vector<std::string> self_trigger_files;
    while((ent = readdir(dir)) != nullptr) {

        std::string name = ent->d_name;
        if(name == "." || name == "..") continue;

        // Consider only .root files
        if(name.size() < 6) continue;
        if(name.substr(name.size()-5) != ".root") continue;

        // Filter external and self-trigger files (both variants)
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if(lower.find("selftrigger") != std::string::npos || lower.find("self_trigger") != std::string::npos) {
            self_trigger_files.push_back(input_dir + name);
            continue;
        }
        external_trigger_files.push_back(input_dir + name);
    }
    closedir(dir);

    // If no files found, print message and exit
    if(external_trigger_files.empty() && self_trigger_files.empty()) {
        std::cout << "No matching ROOT files found in " << input_dir << std::endl;
        return;
    }

    // ------------------------------------------------------------------------------------------------------------
    // Extract HV values and scanned layer information from filenames (assuming format like LY1HV4600)
    std::regex pattern_lyhv("LY(\\d+)HV(\\d+)", std::regex::icase);

    std::map<int, std::set<int>> layer_hvs;                 // Map from layer number to set of HVs encountered (to detect which layer is scanned)
    std::vector<std::pair<std::string, int>> file_and_hv;   // For quick per-file HV lookup for chosen layer

    for(const auto &file_path: external_trigger_files) {

        std::string file_name = file_path.substr(file_path.find_last_of('/') + 1);
        std::smatch match;
        auto begin = file_name.cbegin();
        while(std::regex_search(begin, file_name.cend(), match, pattern_lyhv)) {
            int layer = std::stoi(match[1].str());
            int hv = std::stoi(match[2].str());
            layer_hvs[layer].insert(hv);
            file_and_hv.push_back({file_path, hv});
            begin = match.suffix().first;
        }
    }

    if(layer_hvs.empty()) {
        std::cerr << "No LY..HV patterns found in filenames; aborting." << std::endl;
        return;
    }

    // ------------------------------------------------------------------------------------------------------------
    // Decide which layer was scanned during data taking = the one with largest number of distinct HV values
    int scanned_layer = -1;
    size_t maxcount = 0;
    for(auto &p: layer_hvs) {
        if(p.second.size() > maxcount) {
            maxcount = p.second.size();
            scanned_layer = p.first;
        }
    }

    if(scanned_layer < 0) {
        std::cerr << "Could not determine scanned layer." << std::endl;
        return;
    }

    // Map for layer numbering convention (Facing the outside vs. facing the Interaction Point)
    std::map<int, int> layer_labels = {
        {3, 0}, // Closest to the interaction point (IP)
        {2, 1},
        {1, 2}  // Farthest from the IP, facing the outside
    };
    // Match layer labling conventions
    int hist_layer_index = layer_labels[scanned_layer];

    std::cout << "Detected scanned layer: LY" << scanned_layer << " (" << maxcount << " distinct HVs)" << std::endl;
    std::cout << "For files up to 28.02.2026 convention is LY1/LY2/LY3 (facing the outside) for layers LY2/LY1/LY0 (facing the IP) respectively." << std::endl;

    // ------------------------------------------------------------------------------------------------------------
    // Now for each file, extract HV of the scanned_layer and extract data from the branches/histograms
    std::regex specific_re("LY" + std::to_string(scanned_layer) + "HV(\\d+)", std::regex::icase);

    std::vector<double> hv_values;
    std::vector<std::string> hv_file_paths;

    std::vector<double> eff_eta1_values;
    std::vector<double> eff_eta1_errors;
    std::vector<double> eff_eta2_values;
    std::vector<double> eff_eta2_errors;

    std::vector<double> cluster_size_eta1_means;
    std::vector<double> cluster_size_eta2_means;
    std::vector<double> cluster_size_eta1_errors;
    std::vector<double> cluster_size_eta2_errors;

    for(const auto &file_path: external_trigger_files) {
        std::string file_name = file_path.substr(file_path.find_last_of('/')+1);
        std::smatch m;
        int hv = -1;
        if(std::regex_search(file_name, m, specific_re)) {
            hv = std::stoi(m[1].str());
        } else {
            // Skip files that don't mention the scanned layer
            continue;
        }

        TFile *root_file = TFile::Open(file_path.c_str(), "READ");
        if(!root_file || root_file->IsZombie()) {
            std::cerr << "Failed to open " << file_path << std::endl;
            continue;
        }

        // Efficiency histograms
        std::string dirpath = "efficiencies_histograms/external_plus_rpc_trigger";
        TDirectory *efficiencies_dir = (TDirectory*)root_file->Get(dirpath.c_str());
        if(!efficiencies_dir) {
            std::cerr << "Directory not found in " << file_path << ": " << dirpath << std::endl;
            root_file->Close();
            continue;
        }

        // The stored histogram is named "eff_eta1_rpc_layerN" but contains four bins
        // in order: eta1, eta2, OR, AND. Read bin1->eta1 and bin2->eta2 (Naming is a BUG in the process_raw_hits).
        std::string hist_name = "eff_eta1_rpc_layer" + std::to_string(hist_layer_index);

        TH1 *h = (TH1*)efficiencies_dir->Get(hist_name.c_str());
        if(!h) {
            std::cerr << "Missing histogram in " << file_path << ": " << hist_name << std::endl;
            root_file->Close();
            continue;
        }

        int nb = h->GetNbinsX();
        if(nb < 2) {
            std::cerr << "Histogram " << hist_name << " has fewer than 2 bins; skipping " << file_path << std::endl;
            root_file->Close();
            continue;
        }

        double val_eta1 = h->GetBinContent(1);
        double val_eta2 = h->GetBinContent(2);

        // Calculate binomial errors with n_samples = 1000
        const double n_samples = 1000.0;
        double err_eta1 = sqrt(val_eta1 * (1.0 - val_eta1) / n_samples);
        double err_eta2 = sqrt(val_eta2 * (1.0 - val_eta2) / n_samples);

        // Extract cluster size means and errors from the clusterization tree
        double mean_cs1 = 0.0, error_cs1 = 0.0;
        double mean_cs2 = 0.0, error_cs2 = 0.0;
        
        TTree *clusterization_tree = (TTree*)root_file->Get("clusterization_data");
        if(clusterization_tree) {
            // For cluster_size_eta1
            std::string tmp_histo_name = std::string("__tmp_cs1_") + std::to_string(hv_values.size());
            clusterization_tree->Draw((std::string("cluster_size_eta1>>") + tmp_histo_name).c_str(), "", "goff");
            TH1 *histo_tmp = (TH1*)gDirectory->Get(tmp_histo_name.c_str());
            if(histo_tmp) {
                mean_cs1 = histo_tmp->GetMean();
                double stddev_cs1 = histo_tmp->GetStdDev();
                int n_entries = (int)histo_tmp->GetEntries();
                
                if(n_entries > 0) {
                    error_cs1 = stddev_cs1 / sqrt(n_entries);
                }
                gDirectory->Delete(tmp_histo_name.c_str());
            }

            // For cluster_size_eta2
            tmp_histo_name = std::string("__tmp_cs2_") + std::to_string(hv_values.size());
            clusterization_tree->Draw((std::string("cluster_size_eta2>>") + tmp_histo_name).c_str(), "", "goff");
            histo_tmp = (TH1*)gDirectory->Get(tmp_histo_name.c_str());
            if(histo_tmp) {
                mean_cs2 = histo_tmp->GetMean();
                double stddev_cs2 = histo_tmp->GetStdDev();
                int n_entries = (int)histo_tmp->GetEntries();
                
                if(n_entries > 0) {
                    error_cs2 = stddev_cs2 / sqrt(n_entries);
                }
                gDirectory->Delete(tmp_histo_name.c_str());
            }
        }

        hv_values.push_back(hv);
        hv_file_paths.push_back(file_path);

        eff_eta1_values.push_back(val_eta1);
        eff_eta1_errors.push_back(err_eta1);
        eff_eta2_values.push_back(val_eta2);
        eff_eta2_errors.push_back(err_eta2);

        cluster_size_eta1_means.push_back(mean_cs1);
        cluster_size_eta2_means.push_back(mean_cs2);
        cluster_size_eta1_errors.push_back(error_cs1);
        cluster_size_eta2_errors.push_back(error_cs2);

        root_file->Close();
    }

    if(hv_values.empty()) {
        std::cerr << "No usable files for scanned layer found." << std::endl;
        return;
    }

    // Sort by HV
    std::vector<int> idx(hv_values.size());
    for(size_t i = 0; i < idx.size(); ++i) idx[i]=i;
    std::sort(idx.begin(), idx.end(), [&](int a,int b) { return hv_values[a] < hv_values[b]; });

    std::vector<double> hv_sorted, eff_eta1_values_sorted, eff_eta1_errors_sorted, eff_eta2_values_sorted, eff_eta2_errors_sorted;
    std::vector<std::string> hv_files_sorted;
    std::vector<double> cluster_size_eta1_means_sorted, cluster_size_eta2_means_sorted;
    for(int i: idx) {
        hv_sorted.push_back(hv_values[i]);
        hv_files_sorted.push_back(hv_file_paths[i]);

        eff_eta1_values_sorted.push_back(eff_eta1_values[i]);
        eff_eta1_errors_sorted.push_back(eff_eta1_errors[i]);
        eff_eta2_values_sorted.push_back(eff_eta2_values[i]);
        eff_eta2_errors_sorted.push_back(eff_eta2_errors[i]);

        cluster_size_eta1_means_sorted.push_back(cluster_size_eta1_means[i]);
        cluster_size_eta2_means_sorted.push_back(cluster_size_eta2_means[i]);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Estimate working point (WP) for eta1 and eta2
    // Procedure: normalize each efficiency curve to its maximum (100%),
    // find the point closest to 90% and take WP = HV_at_point + 300 V.
    // Return pair {WP_value, index_in_sorted_arrays}
    auto estimate_wp = [&](const std::vector<double>& eff)->std::pair<double,int> {

        if(eff.empty() || hv_sorted.empty()) return {-1.0, -1};

        double max_eff = *std::max_element(eff.begin(), eff.end());
        if(max_eff <= 0) return {-1.0, -1};

        double target_frac = 0.9; // 90% of normalized curve
        // Find HV corresponding to 90% point (index based on existing HV grid)
        size_t idx90 = 0;
        double best_diff = std::fabs(eff[0] / max_eff - target_frac);
        for(size_t i = 1; i < eff.size(); ++i) {
            double diff = std::fabs(eff[i] / max_eff - target_frac);
            if(diff < best_diff) { best_diff = diff; idx90 = i; }
        }

        double hv_at_90 = hv_sorted[idx90];
        double wp_value = hv_at_90 + 300.0; // nominal working point

        // Choose the file whose HV is closest to wp_value
        int best_idx = 0;
        double best_hvdiff = std::fabs(hv_sorted[0] - wp_value);
        for(size_t j = 1; j < hv_sorted.size(); ++j) {
            double d = std::fabs(hv_sorted[j] - wp_value);
            if(d < best_hvdiff) { best_hvdiff = d; best_idx = (int)j; }
        }

        return { wp_value, best_idx };
    };

    auto res1 = estimate_wp(eff_eta1_values_sorted);
    auto res2 = estimate_wp(eff_eta2_values_sorted);
    double eta1_wp = res1.first; int eta1_idx = res1.second;
    double eta2_wp = res2.first; int eta2_idx = res2.second;

    if(eta1_wp > 0) std::cout << "Working point (eta1) ~ " << eta1_wp << " V (closest HV " << (eta1_idx >= 0 ? hv_sorted[eta1_idx] : -1) << ")" << std::endl;
    else std::cout << "Working point (eta1) could not be determined" << std::endl;
    if(eta2_wp > 0) std::cout << "Working point (eta2) ~ " << eta2_wp << " V (closest HV " << (eta2_idx >= 0 ? hv_sorted[eta2_idx] : -1) << ")" << std::endl;
    else std::cout << "Working point (eta2) could not be determined" << std::endl;


    // ============================================================================================================
    // Plotting session
    gROOT->LoadMacro("AtlasUtils.C");
    gROOT->SetStyle("ATLAS");

    // ------------------------------------------------------------------------------------------------------------
    // Efficiency Curves
    int n_data_points = (int)hv_sorted.size();
    std::vector<double> xerr(n_data_points, 0.0);

    TGraphErrors *g1 = new TGraphErrors(n_data_points, hv_sorted.data(), eff_eta1_values_sorted.data(), xerr.data(), eff_eta1_errors_sorted.data());
    TGraphErrors *g2 = new TGraphErrors(n_data_points, hv_sorted.data(), eff_eta2_values_sorted.data(), xerr.data(), eff_eta2_errors_sorted.data());

    TCanvas *c = new TCanvas("c_eff", "Efficiency vs. High Voltage", 800, 600);

    g1->SetTitle("Efficiency vs. High Voltage");
    g1->GetXaxis()->SetTitle("HV [V]");
    g1->GetYaxis()->SetTitle("Efficiency");
    g1->Draw("APZ");

    // Create a text object for "ATLAS"
    TLatex *atlas = new TLatex();
    atlas->SetNDC();                        // Use normalized coordinates (0-1)
    atlas->SetTextFont(72);                 // Bold font
    atlas->SetTextSize(0.05);               // Adjust size as needed
    atlas->SetTextAlign(12);                // Left-aligned
    atlas->DrawLatex(0.2, 0.85, "ATLAS");

    // Create a text object for "Work in Progress"
    TLatex *work = new TLatex();
    work->SetNDC();
    work->SetTextFont(42);                  // Regular font
    work->SetTextSize(0.04);
    work->SetTextAlign(12);
    work->DrawLatex(0.375, 0.85, "Work in Progress");

    TLine *line = new TLine(5700, gPad->GetUymin(), 5700, gPad->GetUymax());
    line->SetLineColor(kGray+2);
    line->SetLineStyle(2);  // 2 = dashed
    line->SetLineWidth(2);
    line->Draw();

    // Marker
    g1->SetMarkerStyle(21); g1->SetMarkerColor(kRed); g1->SetMarkerSize(1.75); g1->SetLineColor(kRed); g1->SetLineWidth(2);
    g2->SetMarkerStyle(8); g2->SetMarkerColor(kBlue); g2->SetMarkerSize(1.75); g2->SetLineColor(kBlue); g2->SetLineWidth(2);
    g2->Draw("PZ same");

    // Margins
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);

    // Labels
    g1->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
    g1->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more
    g2->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
    g2->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more

    // Ranges
    g1->GetYaxis()->SetRangeUser(0.00, 1.00);
    g2->GetYaxis()->SetRangeUser(0.00, 1.00);

    // Legend
    TLegend *legend = new TLegend(0.65, 0.2, 0.88, 0.35);
    legend->AddEntry(g1, "#eta-", "p");
    legend->AddEntry(g2, "#eta+", "p");
    legend->Draw();

    c->Update();


    // ------------------------------------------------------------------------------------------------------------
    // Cluster mean vs. HV curves

    // Create temporary vectors for HV > 4900
    std::vector<double> hv_filtered;
    std::vector<double> cs1_filtered;
    std::vector<double> cs2_filtered;

    for(int idx = 0; idx < n_data_points; idx++) {
        if(hv_sorted[idx] > 4900) {
            hv_filtered.push_back(hv_sorted[idx]);
            cs1_filtered.push_back(cluster_size_eta1_means_sorted[idx]);
            cs2_filtered.push_back(cluster_size_eta2_means_sorted[idx]);
        }
    }

    // Create graphs with filtered data
    TGraph *cg1 = new TGraphErrors(hv_filtered.size(), hv_filtered.data(), cs1_filtered.data(), nullptr, cluster_size_eta1_errors.data());
    TGraph *cg2 = new TGraphErrors(hv_filtered.size(), hv_filtered.data(), cs2_filtered.data(), nullptr, cluster_size_eta2_errors.data());

    TCanvas *cluster_means_canvas = new TCanvas("c_cluster_means", "Cluster size means vs. High Voltage", 1200, 800);

    cg1->SetTitle("Cluster size mean vs. High voltage");
    cg1->GetXaxis()->SetTitle("HV [V]");
    cg1->GetYaxis()->SetTitle("Mean cluster size");
    cg1->Draw("AP");
    cg2->Draw("P same");

    // Create a text object for "ATLAS"
    atlas->SetNDC();                        // Use normalized coordinates (0-1)
    atlas->SetTextFont(72);                 // Bold font
    atlas->SetTextSize(0.05);               // Adjust size as needed
    atlas->SetTextAlign(12);                // Left-aligned
    atlas->DrawLatex(0.2, 0.85, "ATLAS");

    // Create a text object for "Work in Progress"
    work->SetNDC();
    work->SetTextFont(42);                  // Regular font
    work->SetTextSize(0.04);
    work->SetTextAlign(12);
    work->DrawLatex(0.375, 0.85, "Work in Progress");

    // Marker
    cg1->SetMarkerStyle(21); cg1->SetMarkerColor(kRed); cg1->SetMarkerSize(1.5); cg1->SetLineColor(kRed); cg1->SetLineWidth(1.5);
    cg2->SetMarkerStyle(8); cg2->SetMarkerColor(kBlue); cg2->SetMarkerSize(1.5); cg2->SetLineColor(kBlue); cg2->SetLineWidth(1.5);
    gStyle->SetEndErrorSize(6);

    // Labels
    cg1->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
    cg1->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more
    cg2->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
    cg2->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more

    // Ranges
    cg1->GetYaxis()->SetRangeUser(1.20, 1.50);
    cg2->GetYaxis()->SetRangeUser(1.20, 1.50);
    cg1->GetXaxis()->SetRangeUser(4900, hv_sorted.back() + 100);
    cg2->GetXaxis()->SetRangeUser(4900, hv_sorted.back() + 100);

    // Legend
    TLegend *legend2 = new TLegend(0.70, 0.20, 0.90, 0.35);
    legend2->AddEntry(cg1, "#eta-", "p");
    legend2->AddEntry(cg2, "#eta+", "p");
    legend2->Draw();

    cluster_means_canvas->Update();

    // ------------------------------------------------------------------------------------------------------------
    // Create cluster size @ WP histograms directly
    TH1F *cs_hist1 = nullptr;
    if(eta1_idx >= 0) {
        TFile *f1 = TFile::Open(hv_files_sorted[eta1_idx].c_str(), "READ");
        if(f1 && !f1->IsZombie()) {
            TTree *t1 = (TTree*)f1->Get("clusterization_data");
            if(t1) {
                std::string tmp_histo_name = "__tmp_cs1_hist";
                t1->Draw((std::string("cluster_size_eta1>>") + tmp_histo_name + "(20,0,20)").c_str(), "", "goff");
                cs_hist1 = (TH1F*)gDirectory->Get(tmp_histo_name.c_str());
                cs_hist1->SetDirectory(0);  // Detach from file
            }
            f1->Close();
        }
    }

    TH1F *cs_hist2 = nullptr;
    if(eta2_idx >= 0) {
        TFile *f2 = TFile::Open(hv_files_sorted[eta2_idx].c_str(), "READ");
        if(f2 && !f2->IsZombie()) {
            TTree *t2 = (TTree*)f2->Get("clusterization_data");
            if(t2) {
                std::string tmp_histo_name = "__tmp_cs2_hist";
                t2->Draw((std::string("cluster_size_eta2>>") + tmp_histo_name + "(20,0,20)").c_str(), "", "goff");
                cs_hist2 = (TH1F*)gDirectory->Get(tmp_histo_name.c_str());
                cs_hist2->SetDirectory(0);  // Detach from file
            }
            f2->Close();
        }
    }

    // Create canvas and draw histograms
    TCanvas *cs_hist_canvas = new TCanvas("c_cs_hist", "Cluster size histogram at WP", 800, 600);

    if(cs_hist1) {
        cs_hist1->SetTitle("Cluster size at Working Point");
        cs_hist1->GetXaxis()->SetTitle("Cluster size");
        cs_hist1->GetYaxis()->SetTitle("Number of Clusters");
        
        // Line and marker
        cs_hist1->SetLineColor(kRed); cs_hist1->SetLineWidth(2);
        cs_hist1->SetMarkerStyle(21); cs_hist1->SetMarkerColor(kRed); cs_hist1->SetMarkerSize(1.75);

        // Ticks
        gPad->SetTickx(1);
        gPad->SetTicky(1);

        // Labels
        cs_hist1->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
        cs_hist1->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more
        cs_hist2->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
        cs_hist2->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more

        // Margins
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.05);

        // Ranges
        cs_hist1->GetXaxis()->SetRangeUser(0, 10);
        cs_hist1->GetYaxis()->SetRangeUser(0, 2700);
        
        cs_hist1->Draw("HIST");
        cs_hist1->Draw("P SAME");
    }

    // Create a text object for "ATLAS"
    atlas->SetNDC();                        // Use normalized coordinates (0-1)
    atlas->SetTextFont(72);                 // Bold font
    atlas->SetTextSize(0.05);               // Adjust size as needed
    atlas->SetTextAlign(12);                // Left-aligned
    atlas->DrawLatex(0.2, 0.85, "ATLAS");

    // Create a text object for "Work in Progress"
    work->SetNDC();
    work->SetTextFont(42);                  // Regular font
    work->SetTextSize(0.04);
    work->SetTextAlign(12);
    work->DrawLatex(0.375, 0.85, "Work in Progress");

    if(cs_hist2) {

        // Line and marker
        cs_hist2->SetLineColor(kBlue); cs_hist2->SetLineWidth(2);
        cs_hist2->SetMarkerStyle(20); cs_hist2->SetMarkerColor(kBlue); cs_hist2->SetMarkerSize(1.5);
        
        cs_hist2->Draw("HIST SAME");
        cs_hist2->Draw("P SAME");
    }

    // Add legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.88, 0.88);
    if(cs_hist1) leg->AddEntry(cs_hist1, "#eta-", "lp");
    if(cs_hist2) leg->AddEntry(cs_hist2, "#eta+", "lp");
    leg->Draw();

    cs_hist_canvas->Update();

    // ------------------------------------------------------------------------------------------------------------
    // Create ToT @ WP histograms directly
    TH1F *tot_hist1 = nullptr;
    if(eta1_idx >= 0) {
        TFile *f1 = TFile::Open(hv_files_sorted[eta1_idx].c_str(), "READ");
        if(f1 && !f1->IsZombie()) {
            TTree *t1 = (TTree*)f1->Get("clusterization_data");
            if(t1) {
                // First create histogram with original binning (40 bins from 0-40 ticks)
                std::string tmp_histo_name = "__tmp_tot1_hist";
                t1->Draw((std::string("cluster_tot1>>") + tmp_histo_name + "(40,0,40)").c_str(), "cluster_tot1>0", "goff");
                TH1F *hist_tmp = (TH1F*)gDirectory->Get(tmp_histo_name.c_str());
                
                if(hist_tmp) {
                    // Create new histogram with nanosecond binning
                    // Convert range: 0-40 ticks → 0-33.32 ns (40 * 0.833 = 33.32)
                    tot_hist1 = new TH1F("tot_hist1", "Time over Threshold;ToT [ns];Entries", 
                                        40, 0, 33.32);
                    
                    // Fill new histogram with converted values
                    for(int i = 1; i <= hist_tmp->GetNbinsX(); i++) {
                        double bin_center_ticks = hist_tmp->GetBinCenter(i);
                        double bin_center_ns = bin_center_ticks * 0.833;
                        double bin_content = hist_tmp->GetBinContent(i);
                        
                        // Find corresponding bin in new histogram and add content
                        int new_bin = tot_hist1->FindBin(bin_center_ns);
                        tot_hist1->AddBinContent(new_bin, bin_content);
                    }
                    
                    // Copy statistical information
                    tot_hist1->SetEntries(hist_tmp->GetEntries());
                    tot_hist1->SetDirectory(0);
                    
                    gDirectory->Delete(tmp_histo_name.c_str());
                }
            }
            f1->Close();
        }
    }

    TH1F *tot_hist2 = nullptr;
    if(eta2_idx >= 0) {
        TFile *f2 = TFile::Open(hv_files_sorted[eta2_idx].c_str(), "READ");
        if(f2 && !f2->IsZombie()) {
            TTree *t2 = (TTree*)f2->Get("clusterization_data");
            if(t2) {
                std::string tmp_histo_name = "__tmp_tot2_hist";
                t2->Draw((std::string("cluster_tot2>>") + tmp_histo_name + "(40,0,40)").c_str(), "cluster_tot2>0", "goff");
                TH1F *hist_tmp = (TH1F*)gDirectory->Get(tmp_histo_name.c_str());
                
                if(hist_tmp) {
                    tot_hist2 = new TH1F("tot_hist2", "Time over Threshold;ToT [ns];Entries", 
                                        40, 0, 33.32);
                    
                    for(int i = 1; i <= hist_tmp->GetNbinsX(); i++) {
                        double bin_center_ticks = hist_tmp->GetBinCenter(i);
                        double bin_center_ns = bin_center_ticks * 0.833;
                        double bin_content = hist_tmp->GetBinContent(i);
                        
                        int new_bin = tot_hist2->FindBin(bin_center_ns);
                        tot_hist2->AddBinContent(new_bin, bin_content);
                    }
                    
                    tot_hist2->SetEntries(hist_tmp->GetEntries());
                    tot_hist2->SetDirectory(0);
                    
                    gDirectory->Delete(tmp_histo_name.c_str());
                }
                f2->Close();
            }
        }
    }

    // Create canvas for ToT histograms
    TCanvas *tot_hist_canvas = new TCanvas("tot_hist_canvas", "ToT Histograms", 800, 600);

    if(tot_hist1) {
        tot_hist1->SetTitle("Time over Threshold");
        tot_hist1->GetXaxis()->SetTitle("Time over Threshold [ns]");
        tot_hist1->GetYaxis()->SetTitle("Number of Hits");
        
        // Line
        tot_hist1->SetLineColor(kRed); tot_hist1->SetLineWidth(2);

        // Ticks
        gPad->SetTickx(1);
        gPad->SetTicky(1);

        // Labels
        tot_hist1->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
        tot_hist1->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more
        tot_hist2->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
        tot_hist2->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more

        // Margins
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.05);

        // Ranges
        tot_hist1->GetXaxis()->SetRangeUser(0, 40);
        tot_hist1->GetYaxis()->SetRangeUser(0, 650);
        
        tot_hist1->Draw("HIST");
    }

    if(tot_hist2) {

        // Line
        tot_hist2->SetLineColor(kBlue); tot_hist2->SetLineWidth(2);
        
        tot_hist2->Draw("HIST SAME");
    }

    // ATLAS WIP label
    atlas->DrawLatex(0.2, 0.85, "ATLAS");
    work->DrawLatex(0.375, 0.85, "Work in Progress");

    // Add legend
    TLegend *leg_tot = new TLegend(0.7, 0.7, 0.88, 0.88);
    if(tot_hist1) leg_tot->AddEntry(tot_hist1, "#eta-", "l");
    if(tot_hist2) leg_tot->AddEntry(tot_hist2, "#eta+", "l");
    leg_tot->Draw();

    tot_hist_canvas->Update();

    // ------------------------------------------------------------------------------------------------------------
    // Create graph and canvas for efficiency vs. v_th

    std::vector<double> vth_values_low_hv = {2.2, 2.3, 2.4, 2.5};
    std::vector<double> eff_vth_low_hv = {0.7686, 0.75, 0.7784, 0.7771};
    std::vector<double> eff_vth_high_hv = {0.9724, 0.96, 0.9733, 0.9567};

    TGraph *g_vth_low = new TGraph(vth_values_low_hv.size(), vth_values_low_hv.data(), eff_vth_low_hv.data());
    TGraph *g_vth_high = new TGraph(vth_values_low_hv.size(), vth_values_low_hv.data(), eff_vth_high_hv.data());

    TCanvas *c_vth = new TCanvas("c_vth", "Efficiency vs. Threshold voltage", 800, 600);

    g_vth_low->SetTitle("Efficiency vs. Threshold voltage");
    g_vth_low->GetXaxis()->SetTitle("Threshold voltage [V]");
    g_vth_low->GetYaxis()->SetTitle("Efficiency");

    g_vth_low->SetMarkerStyle(21); g_vth_low->SetMarkerColor(kRed); g_vth_low->SetMarkerSize(1.75); g_vth_low->SetLineColor(kRed); g_vth_low->SetLineWidth(2);
    g_vth_high->SetMarkerStyle(8); g_vth_high->SetMarkerColor(kBlue); g_vth_high->SetMarkerSize(1.75); g_vth_high->SetLineColor(kBlue); g_vth_high->SetLineWidth(2);

    g_vth_low->Draw("AP");
    g_vth_high->Draw("P SAME");

    atlas->DrawLatex(0.2, 0.85, "ATLAS");
    work->DrawLatex(0.375, 0.85, "Work in Progress");

    // Margins
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);

    g_vth_low->GetXaxis()->SetRangeUser(2.15, 2.55);
    g_vth_low->GetYaxis()->SetRangeUser(0.7, 1.0);
    g_vth_high->GetXaxis()->SetRangeUser(2.15, 2.55);
    g_vth_high->GetYaxis()->SetRangeUser(0.7, 1.0);

    // Legend
    TLegend *legend_vth = new TLegend(0.65, 0.2, 0.88, 0.35);
    legend_vth->AddEntry(g_vth_low, "5.3 kV", "p");
    legend_vth->AddEntry(g_vth_high, "5.8 kV", "p");
    legend_vth->Draw();

    // ------------------------------------------------------------------------------------------------------------
    // Create graph and canvas for noise rate vs. HV (from Giorgia's data)

    std::vector<double> noise_hvs_giorgia = {4600, 4800, 5000, 5200, 5800};
    std::vector<double> noise_rates_giorgia = {0.0001, 0.1220, 0.2470, 0.3730, 0.5000};
    std::vector<double> efficiencies = {0.001, 0.016, 0.229, 0.712, 0.956};
    std::vector<int> n_hits = {};

    int n_points = noise_hvs_giorgia.size();
    
    // Arrays for TGraphErrors
    std::vector<double> hv_errors(n_points, 0.0); // No uncertainty on the x-axis (HV)
    std::vector<double> noise_errors(n_points, 0.0);

    double time_window = 400e-9;         // 400 ns in seconds
    int n_strips = 48;
    double area_strip = 2.5 * 180.0;     // cm^2
    double total_area = area_strip * n_strips;

    // Calculate Poisson errors
    for (int i = 0; i < n_points; ++i) {
        double n_events = 2000.0 * efficiencies[i];

        // Prevent division by zero if efficiency is 0.000 (e.g., at 4600V)
        if (n_events > 0 && noise_rates_giorgia[i] > 0) {

            // Back-calculate the actual number of hits (N)
            // rate = N / (time_window * total_area * n_events)  --> N = rate * (time_window * total_area * n_events)
            double n_hits = noise_rates_giorgia[i] * (time_window * total_area * n_events);

            // Error on N is sqrt(N)
            double n_hits_err = std::sqrt(n_hits);

            // Error on rate as sqrt(N) / (time_window * total_area * n_events)
            noise_errors[i] = n_hits_err / (time_window * total_area * n_events);

            std::cout << "HV: " << noise_hvs_giorgia[i] << " V, Noise rate: " << noise_rates_giorgia[i] 
                      << " Hz/cm^2, N_hits: " << n_hits << ", Rate error: " << noise_errors[i] << std::endl;

        } else {
            // Can't reliably back-calculate the error with this formula if efficiency is zero
            noise_errors[i] = 0.0;
        }
    }

    // Graph and Canvas
    TGraphErrors *g_noise = new TGraphErrors(n_points, noise_hvs_giorgia.data(), noise_rates_giorgia.data(), hv_errors.data(), noise_errors.data());

    TCanvas *c_noise_giorgia = new TCanvas("c_noise_giorgia", "Noise Rate vs. High Voltage", 800, 600);

    g_noise->SetTitle("Noise rate vs. High Voltage");
    g_noise->GetXaxis()->SetTitle("HV [V]");
    g_noise->GetYaxis()->SetTitle("Noise rate [Hz/cm^{2}]");

    g_noise->SetMarkerStyle(8); g_noise->SetMarkerColor(kBlack); g_noise->SetMarkerSize(1.5); 
    g_noise->SetLineColor(kBlack); g_noise->SetLineWidth(2);

    g_noise->Draw("AP");

    atlas->DrawLatex(0.2, 0.85, "ATLAS");
    work->DrawLatex(0.375, 0.85, "Work in Progress");

    // Margins
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);

    // Legend
    TLegend *legend_g = new TLegend(0.65, 0.2, 0.88, 0.35);
    legend_g->AddEntry(tot_hist1, "#eta-", "pl");
    legend_g->Draw();

    c_noise_giorgia->Update();

    // ------------------------------------------------------------------------------------------------------------
    // Copy clusterization_data at WP into output file

    std::string output_file = "analysis_" + date_str + "_LY" + std::to_string(scanned_layer) + ".root";
    TFile of(output_file.c_str(), "RECREATE");

    // Copy clusterization_data tree from the file closest to WP for eta1 and eta2
    auto copy_tree_from = [&](int idx, const std::string &label) {
        if(idx < 0 || idx >= (int)hv_files_sorted.size()) return;
        std::string srcpath = hv_files_sorted[idx];
        TFile *src = TFile::Open(srcpath.c_str(), "READ");
        if(!src || src->IsZombie()) {
            std::cerr << "Failed to open source file for copying: " << srcpath << std::endl;
            if(src) src->Close();
            return;
        }
        TTree *t = (TTree*)src->Get("clusterization_data");
        if(!t) {
            std::cerr << "No clusterization_data tree in " << srcpath << std::endl;
            src->Close();
            return;
        }
        // Write a cloned tree into the output file under a name indicating which WP it corresponds to
        std::string outname = std::string("clusterization_data_") + label + "_WP";
        of.cd();
        TTree *tclone = t->CloneTree(-1, "fast");
        if(tclone) {
            tclone->SetName(outname.c_str());
            tclone->Write(outname.c_str());
            delete tclone;
            std::cout << "Copied " << outname << " from " << srcpath << std::endl;
        } else {
            std::cerr << "Failed to clone tree from " << srcpath << std::endl;
        }
        src->Close();
    };

    // If both indices point to same file, copy only once
    if(eta1_idx == eta2_idx) {
        if(eta1_idx >= 0) copy_tree_from(eta1_idx, "eta1_eta2");
    } else {
        if(eta1_idx >= 0) copy_tree_from(eta1_idx, "eta1");
        if(eta2_idx >= 0) copy_tree_from(eta2_idx, "eta2");
    }

    // ============================================================================================================
    // Noise rate vs. High voltage

    std::vector<double> noise_hvs;
    std::vector<double> noise_rates;
    std::vector<double> noise_rate_errors;

    for (auto &file: self_trigger_files) {
        std::string file_name = file.substr(file.find_last_of('/')+1);
        std::smatch match2;
        int hv = -1;
        if(std::regex_search(file_name, match2, specific_re)) {
            hv = std::stoi(match2[1].str());
        } else continue;

        TFile *root_file = TFile::Open(file.c_str(), "READ");
        if(!root_file || root_file->IsZombie()) { if(root_file) root_file->Close(); continue; }
        TTree *processed_data_tree = (TTree*)root_file->Get("processed_data");
        if(!processed_data_tree) { root_file->Close(); continue; }

        Long64_t total_events_hv = processed_data_tree->GetEntries();
        if(total_events_hv <= 0) { root_file->Close(); continue; }

        // Count hits in the selected layer across all strips
        std::string selection = std::string("proc_layer==") + std::to_string(hist_layer_index);
        std::string tmpname = std::string("__tmp_noise_strip_hv") + std::to_string(hv);
        processed_data_tree->Draw((std::string("proc_strip>>") + tmpname).c_str(), selection.c_str(), "goff");
        TH1 *histo_tmp = (TH1*)gDirectory->Get(tmpname.c_str());

        double total_hits = 0.0;
        if(histo_tmp) {
            total_hits = histo_tmp->GetEntries();
            gDirectory->Delete(tmpname.c_str());
        }

        // Rate = hits / (time_window_s * area_strip * n_strips * total_events)
        double rate = total_hits / (time_window * area_strip * (double)n_strips * (double)total_events_hv);

        // Statistical error: assuming Poisson statistics on hits count
        // sigma_rate = sqrt(hits) / (time_window_s * area_strip * n_strips * total_events)
        double rate_error = (total_hits > 0.0)
            ? std::sqrt(total_hits) / (time_window * area_strip * (double)n_strips * (double)total_events_hv)
            : 0.0;

        noise_hvs.push_back((double)hv);
        noise_rates.push_back(rate);
        noise_rate_errors.push_back(rate_error);

        std::cout << "HV " << hv << " V: " << total_hits << " hits in " << total_events_hv
                  << " events => noise rate = " << rate << " +/- " << rate_error << " Hz/cm^2" << std::endl;

        root_file->Close();
    }

    // Sort by HV
    if(!noise_hvs.empty()) {
        std::vector<int> noise_idx(noise_hvs.size());
        for(size_t i = 0; i < noise_idx.size(); ++i) noise_idx[i] = i;
        std::sort(noise_idx.begin(), noise_idx.end(), [&](int a, int b){ return noise_hvs[a] < noise_hvs[b]; });

        std::vector<double> noise_hvs_s, noise_rates_s, noise_rate_errors_s;
        for(int i: noise_idx) {
            noise_hvs_s.push_back(noise_hvs[i]);
            noise_rates_s.push_back(noise_rates[i]);
            noise_rate_errors_s.push_back(noise_rate_errors[i]);
        }
        noise_hvs        = noise_hvs_s;
        noise_rates      = noise_rates_s;
        noise_rate_errors = noise_rate_errors_s;
    }

    // ------------------------------------------------------------------------------------------------------------
    // Create graph and canvas
    TGraphErrors *noise_hv_graph = nullptr;
    TCanvas *noise_hv_canvas = nullptr;

    if(!noise_hvs.empty()) {
        std::vector<double> noise_hv_xerr(noise_hvs.size(), 0.0);
        noise_hv_graph = new TGraphErrors(
            (int)noise_hvs.size(),
            noise_hvs.data(),
            noise_rates.data(),
            noise_hv_xerr.data(),
            noise_rate_errors.data()
        );

        noise_hv_canvas = new TCanvas("c_noise_hv", "Noise rate vs. High Voltage", 800, 600);
        noise_hv_graph->SetTitle("Noise rate vs. High Voltage");
        noise_hv_graph->GetXaxis()->SetTitle("HV [V]");
        noise_hv_graph->GetYaxis()->SetTitle("Noise rate [Hz/cm^{2}]");

        noise_hv_graph->SetMarkerStyle(20); noise_hv_graph->SetMarkerSize(1.5); noise_hv_graph->SetMarkerColor(kRed);
        noise_hv_graph->SetLineColor(kRed); noise_hv_graph->SetLineWidth(2); gStyle->SetEndErrorSize(6);

        noise_hv_graph->Draw("APZ");

        // Margins
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.05);

        noise_hv_graph->GetXaxis()->SetTitleOffset(1.2);
        noise_hv_graph->GetYaxis()->SetTitleOffset(1.4);

        // ATLAS WIP label
        atlas->DrawLatex(0.2, 0.85, "ATLAS");
        work->DrawLatex(0.375, 0.85, "Work in Progress");

        noise_hv_canvas->Update();
    }


    // ============================================================================================================
    // Noise rate per strip calculation at the working point using self-trigger files
    // Choose the self-trigger file HV using ceil behavior: the smallest available self-trigger HV >= nominal WP
    TGraph *noise_graph = nullptr;
    TCanvas *noise_canvas = nullptr;

    int total_events = 0;
    std::vector<int> strip_counts(n_strips, 0);

    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Starting noise rate calculation at (or next to) working point..." << std::endl;

    // Determine target HV for noise calculation based on nominal WP and available self-trigger files
    int target_hv = -1;
    double nominal_wp = (eta1_wp > 0.0) ? eta1_wp : ((eta2_wp > 0.0) ? eta2_wp : -1.0);
    if(nominal_wp <= 0.0) {
        std::cout << "\nNo working point available; skipping noise calculation." << std::endl;
    } else {
        // Get all available HVs from self-trigger files
        std::vector<int> self_hvs;
        for(const auto &sf: self_trigger_files) {
            std::string fn = sf.substr(sf.find_last_of('/')+1);
            std::smatch m2;
            if(std::regex_search(fn, m2, specific_re)) {
                self_hvs.push_back(std::stoi(m2[1].str()));
            }
        }
        if(self_hvs.empty()) {
            std::cout << "\nNo self-trigger files with HV found; skipping noise calculation." << std::endl;
        } else {
            std::sort(self_hvs.begin(), self_hvs.end());
            self_hvs.erase(std::unique(self_hvs.begin(), self_hvs.end()), self_hvs.end());
            int ceil_val = (int)std::ceil(nominal_wp);
            auto it = std::lower_bound(self_hvs.begin(), self_hvs.end(), ceil_val);
            target_hv = (it == self_hvs.end()) ? self_hvs.back() : *it;
            std::cout << "\n";
            std::cout << "Nominal WP: " << nominal_wp << " V; Selected self-trigger HV: " << target_hv << " V (Ceil behavior)" << std::endl;

            for(const auto &file: self_trigger_files) {

                // Find the file corresponding to the target HV for the scanned layer
                std::string file_name = file.substr(file.find_last_of('/')+1);
                std::smatch match2;
                if(std::regex_search(file_name, match2, specific_re)) {
                    int hv = std::stoi(match2[1].str());
                    if(hv != target_hv) continue;
                } else continue;

                TFile *root_file = TFile::Open(file.c_str(), "READ");
                if(!root_file || root_file->IsZombie()) { if(root_file) root_file->Close(); continue; }
                TTree *processed_data_tree = (TTree*)root_file->Get("processed_data");
                if(!processed_data_tree) { root_file->Close(); continue; }

                // Use TTree::Draw to histogram proc_strip for the selected layer
                Long64_t entries = processed_data_tree->GetEntries();
                total_events += (int)entries;
                std::string tmpname = std::string("__tmp_strip_") + std::to_string(total_events);
                std::string selection = std::string("proc_layer==") + std::to_string(hist_layer_index);
                processed_data_tree->Draw((std::string("proc_strip>>") + tmpname).c_str(), selection.c_str(), "goff");
                TH1 *histo_tmp = (TH1*)gDirectory->Get(tmpname.c_str());
                if(histo_tmp) {
                    int nbins = histo_tmp->GetNbinsX();
                    for(int b = 1; b <= nbins; ++b) {
                        double center = histo_tmp->GetXaxis()->GetBinCenter(b);
                        std::cout << "Bin " << b << ": center=" << center;
                        int strip_idx = (int)std::floor(center + 1e-9);
                        double content = histo_tmp->GetBinContent(b);
                        if(strip_idx >= 0 && strip_idx < n_strips) {
                            strip_counts[strip_idx] += (int)content;
                        }
                    }
                    gDirectory->Delete(tmpname.c_str());
                }

                root_file->Close();
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------------
    // Compute rates: hits / (time_window * area_strip * n_events)
    std::vector<double> strip_rates(n_strips, 0.0);
    if(total_events > 0) {
        for(int idx = 0; idx < n_strips; ++idx) {
            double hits = (double)strip_counts[idx];
            std::cout << "Strip " << idx << ": " << hits << " hits in " << total_events << " events." << std::endl;
            strip_rates[idx] = hits / (time_window * area_strip * (double)total_events);
        }
    }

    // ------------------------------------------------------------------------------------------------------------
    // Create graph and canvas
    std::vector<double> strip_x(n_strips);
    for(int idx = 0; idx < n_strips; ++idx) strip_x[idx] = idx;

    TGraph *ng = new TGraph(n_strips, strip_x.data(), strip_rates.data());
    ng->SetMarkerStyle(20); ng->SetMarkerColor(kRed);

    TCanvas *nc = new TCanvas("c_noise", "Noise rate per strip", 800, 600);
    ng->SetTitle((std::string("Noise rate per strip at HV ") + std::to_string(target_hv) + " V - LY" + std::to_string(scanned_layer)).c_str());
    ng->GetXaxis()->SetTitle("Strip");
    ng->GetYaxis()->SetTitle("Noise rate [Hz/cm^{2}]");
    ng->Draw("AP");

    // Marker
    ng->SetMarkerStyle(20);
    ng->SetMarkerSize(1.5);
    ng->SetMarkerColor(kRed);

    // Ticks
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    
    // Grid
    gPad->SetGrid();

    // Margins
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);

    // Labels
    ng->GetXaxis()->SetTitleOffset(1.2);  // Default is 1.0, increase for more space
    ng->GetYaxis()->SetTitleOffset(1.4);  // Y usually needs more

    // Ranges
    ng->GetXaxis()->SetRangeUser(-0.5, 47.5);

    nc->Update();

    noise_graph = ng;
    noise_canvas = nc;

    // Save noise rate graph
    of.cd();
    of.mkdir("noise_rates");
    of.cd("noise_rates");
    noise_graph->Write("noise_rate_graph");

    // ------------------------------------------------------------------------------------------------------------
    // Finally store all the graphs and canvases

    // Save efficiency graphs into output file under efficiency_curves
    of.mkdir("efficiency_curves");
    of.cd("efficiency_curves");
    g1->Write("eff_eta1_graph");
    g2->Write("eff_eta2_graph");

    // Save cluster means graphs into output file under cluster_means
    of.cd();
    of.mkdir("cluster_means");
    of.cd("cluster_means");
    if(cg1) { cg1->Write("cluster_size_eta1_means_graph"); }
    if(cg2) { cg2->Write("cluster_size_eta2_means_graph"); }

    // Finally save all the canvases into output file under canvases
    of.cd();
    of.mkdir("canvases");
    of.cd("canvases");
    c->Write("canvas_eff_overlay");
    cluster_means_canvas->Write("canvas_cluster_means");
    noise_canvas->Write("noise_rate_canvas");
    cs_hist_canvas->Write("canvas_cluster_size_hist");
    tot_hist_canvas->Write("canvas_tot_hist");
    noise_hv_canvas->Write("noise_rate_hv_canvas");
    c_vth->Write("canvas_eff_vth");
    c_noise_giorgia->Write("canvas_noise_giorgia");

    of.Close();

    std::cout << "Saved analysis graphs and cluster data to " << output_file << std::endl;
}