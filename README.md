# ROME1 DCT Efficiency Analysis

Analysis framework for RPC certification/testing with the ROME1-DCT.

## Quick Start

### Build

```bash
make clean
make build
```

### Data Acquisition

```bash
./scripts/run_data_acquisition.sh <nEvents> <comment> [OPTIONS]
```

**Examples:**
```bash
# Basic batch mode acquisition
./scripts/run_data_acquisition.sh 1000 test_run -b

# Self-trigger mode
./scripts/run_data_acquisition.sh 1000 test_run -s --batch
```

**Options:**
- `-b`, `--batch` — Run Vivado in batch mode
- `-s`, `--self` — Use self-trigger firmware
- `-h`, `--help` — Display help message

### Analysis

```bash
./scripts/run_analysis.sh <input> [OPTIONS]
```

**Examples:**
```bash
# Analyze a single raw data file
./scripts/run_analysis.sh data/raw/2026-01-28_11-55.txt --dt-max -120 --dt-min -200

# Batch process all files in a directory
./scripts/run_analysis.sh data/raw/ --self --dt-max -150
```

**Options:**
- `-s`, `--self` — Used self-trigger firmware during acquisition
- `--dt-max VALUE` — Maximum time window for efficiency (default: -100)
- `--dt-min VALUE` — Minimum time window for efficiency (default: -180)
- `-h`, `--help` — Display help message

## Directory Structure

```
├── bin/                          # Compiled executables
│   └── analysis                  # Main executable for hits processing and analysis
├── data/
│   ├── raw/                      # Raw data text files from data acquisition
│   └── output/                   # Generated ROOT files and plots
├── include/
│   ├── cluster.hpp               # Cluster data structure and methods
│   ├── constants.hpp             # Analysis constants and configuration
│   ├── dataAnalyzer.hpp          # Main analyzer class
│   ├── event.hpp                 # Event data structure
│   ├── hit.hpp                   # Hit data structure
│   ├── track.hpp                 # Track reconstruction
│   └── types.hpp                 # Type definitions
├── src/
│   ├── main.cpp                  # Entry point and orchestration
│   ├── dataAnalyzer.cpp          # Core analysis implementation
│   ├── cluster.cpp               # Cluster processing
│   ├── event.cpp                 # Event handling
│   ├── hit.cpp                   # Hit processing
│   └── track.cpp                 # Track reconstruction
├── utils/
│   └── utils.hpp                 # General utility functions
├── scripts/
│   ├── run_data_acquisition.sh   # Data acquisition orchestrator
│   ├── run_analysis.sh           # Batch analysis on raw data files
│   └── GO.tcl                    # Vivado TCL script for data acquisition
└── Makefile                      # Build system (C++17, ROOT 6.x)
```

## Processing Pipeline

**Phase 1: Data Acquisition**
- 1. Vivado acquires binary data (`tmp_file*`)
- 2. Bash merge & filter → Raw data file (`*.txt`)

**Phase 2: Analysis**
- 3. C++ processing on each raw file → Efficiency calculation and histograms
- 4. Output → ROOT trees and PDF plots

## Time Window Parameters

- `dt_max` — Maximum time difference for efficiency window (negative values)
- `dt_min` — Minimum time difference for efficiency window (more negative)
- Efficiency calculated within: `dt_min ≤ dt ≤ dt_max`

## Dependencies

- ROOT 6.x (TFile, TTree, TH1)
- C++17 compiler (clang++ on macOS)
- Vivado (for data acquisition)
- bash (for data merging and filtering)
