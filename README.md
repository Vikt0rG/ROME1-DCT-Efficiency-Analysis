# ROME1 DCT Efficiency Analysis

Analysis framework for RPC certification/testing with the ROME1-DCT.

## Quick Start

### Build

```bash
make clean
make process_raw_hits
```

### Run Data Acquisition & Processing

```bash
./scripts/run_data_acquisition.sh <nEvents> <comment> [OPTIONS]
```

**Examples:**
```bash
# Basic batch mode acquisition
./scripts/run_data_acquisition.sh 1000 test_run -b

# Self-trigger with custom time windows
./scripts/run_data_acquisition.sh 1000 test_run -s --dt-max -120 --dt-min -200

# Shorthand and longform options are equivalent
./scripts/run_data_acquisition.sh 1000 test_run --batch --self --dt-max -150
```

**Options:**
- `-b`, `--batch` — Run Vivado in batch mode
- `-s`, `--self` — Use self-trigger firmware
- `--dt-max VALUE` — Maximum time window for efficiency (default: -100)
- `--dt-min VALUE` — Minimum time window for efficiency (default: -180)
- `-h`, `--help` — Display help message

## Directory Structure

```
├── bin/                          # Compiled executables
│   └── process_raw_hits          # Standalone binary for raw hits processing
├── data/
│   ├── input/                    # Input ROOT files organized by date
│   └── output/                   # Generated ROOT files and plots
├── include/
│   └── process_raw_hits.hpp      # Header with function declarations
├── src/
│   ├── main.cpp                  # Entry point and orchestration
│   ├── process_raw_hits_utils.cpp # Modular processing functions
│   └── utils.cpp                 # General utilities
├── root_macros/                  # ROOT analysis macros
├── scripts/
│   ├── run_data_acquisition.sh   # Main workflow orchestrator
│   ├── GO.tcl                    # Vivado TCL script
│   └── merge_ila_files.py        # Python data merge utility
└── Makefile                      # Build system (C++17, ROOT 6.x)
```

## Processing Pipeline

1. **Vivado Acquisition** → Binary data (`tmp_file*`)
2. **Python Merge** → ILA data file (`iladata.txt`)
3. **Strip & Filter** → Processed data (`tmp.strip`)
4. **C++ Processing** → Efficiency calculation and histograms
5. **Output** → ROOT trees and PDF plots

## Time Window Parameters

- `dt_max` — Maximum time difference for efficiency window (negative values)
- `dt_min` — Minimum time difference for efficiency window (more negative)
- Efficiency calculated within: `dt_min ≤ dt ≤ dt_max`

## Dependencies

- ROOT 6.x (TFile, TTree, TH1)
- C++17 compiler (clang++ on macOS)
- Vivado (for data acquisition)
- Python 3 (for data processing)
