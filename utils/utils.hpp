#pragma once

#include <string>
#include <vector>
#include <TDirectory.h>
#include <TFile.h>

/// @namespace TimeUtils
/// @brief Namespace for utility functions related to time conversions and calculations
namespace TimeUtils {
[[maybe_unused]] static int convertRawTimeToPhysical(int raw_time, int bcid) {
    return raw_time != 0 ? bcid * 30 + raw_time - 1 : -1;
}
}

/// @namespace PathUtils
/// @brief Namespace for utility functions related to path (and directory) management
namespace PathUtils {

/// @brief Utility function to ensure a directory exists in the ROOT file and return a
/// pointer to it
/// @param parent Pointer to the parent directory where the new directory should be created
/// @param name Pointer to the name of the directory to be created
/// @return Pointer to the created or existing directory, or nullptr if there was an error 
TDirectory* ensureDirectory(TDirectory* parent, const char* name) {
    TDirectory* dir = parent->GetDirectory(name);
    if (!dir) {
        dir = parent->mkdir(name);
    }
    return dir;
}

/// @brief Overloaded function to make sure a directory exists in the ROOT file and return
/// a pointer to it
/// @param file Pointer to the ROOT file where the new directory should be created
/// @param name Pointer to the name of the directory to be created
/// @return Pointer to the created or existing directory, or nullptr if there was an error
TDirectory* ensureDirectory(TFile* file, const char* name) {
    TDirectory* dir = file->GetDirectory(name);
    if (!dir) {
        dir = file->mkdir(name);
    }
    return dir;
}

/// @brief Utility function to setup the directory structure in the input ROOT file
/// @param input_file Pointer to the input ROOT file where the directory structure will be
/// created
/// @param dir_names Vector of directory paths to be created (e.g. "analysis/ToT_vs_Strip")
void setupDirectories(TFile* input_file, const std::vector<std::string>& dir_names) {
    for (const auto& full_path : dir_names) {
        std::stringstream ss(full_path);
        std::string component;
        
        // Start at the root file level
        TDirectory* current_dir = input_file; 
        
        // Split using "/" as a delimiter and create/access directories iteratively
        while (std::getline(ss, component, '/')) {
            if (component.empty()) continue; // Skip accidental double slashes "//"
            
            // Use TDirectory* recursively to ensure each level of the path exists
            current_dir = ensureDirectory(current_dir, component.c_str());
        }
    }
}

}