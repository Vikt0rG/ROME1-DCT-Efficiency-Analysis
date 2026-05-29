#pragma once

#include <map>
#include <string>
#include <vector>

#include <TTree.h>
#include <TBranch.h>
#include <TF1.h>
#include <TH1.h>
#include <TFitResult.h>

/// @brief Reads data from a TTree branch and fills a histogram.
/// @param tree 
/// @param branch 
/// @return Pointer to the filled histogram, or nullptr on error.
TH1* readDataAsHist(TTree* tree, TBranch* branch);

/// @brief Creates a TF1 representing a flat background plus a Gaussian signal.
/// @return Pointer to the created TF1, or nullptr on error.
TF1* flatPlusGaussian();

/// @brief Fits the provided histogram with the given fit function and initial parameters.
/// @param hist 
/// @param fitFunc 
/// @param initialParams 
/// @return Vector of fitted parameters, or empty vector on error.
std::vector<double> fitBackground(TH1* hist, TF1* fitFunc, const std::vector<double> initialParams);

/// @brief Sets the rejection range based on the fitted function parameters.
/// @param fitFunc 
/// @return Pair of integers representing the signal range.
std::pair<int, int> setSignalRange(TF1* fitFunc);

/// @struct SelectionMask
/// @brief Represents a selection mask with a name and a boolean vector indicating selected entries.
/// @param name Name of the selection mask.
/// @param mask Boolean vector where true indicates selected entries.
struct SelectionMask;

/// @brief Generates selection masks for rising and falling edges based on hit times and rise information.
/// @param hit_time1 
/// @param hit_time2 
/// @param hit_rise
/// @return Vector of SelectionMask objects for different selections.
std::vector<SelectionMask> getMasks(const std::vector<int>& hit_time1, const std::vector<int>& hit_time2, const std::vector<int>& hit_rise);

/// @brief Creates a new branch in the TTree to flag signal events based on defined category ranges.
/// @param tree Pointer to the TTree to modify.
/// @param categoryRanges Map of category names to their corresponding low and high range values for selection
void createSignalBranch(TTree* tree, const std::map<std::string, std::pair<double, double>>& categoryRanges);