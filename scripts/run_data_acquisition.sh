#!/bin/bash

# Script to run the full data acquisition and processing workflow for the DCT test setup

# Check for help flag first
if [[ "$@" == *"--help"* ]]; then
   echo "Script usage: $0 <nEvents> <comment> [OPTIONS]"
   echo ""
   echo "REQUIRED ARGUMENTS:"
   echo "  <nEvents>           Number of events to acquire"
   echo "  <comment>           Comment/description for the run"
   echo ""
   echo "OPTIONS:"
   echo "  -b, --batch         Run Vivado in batch mode"
   echo "  -s, --self          Use self-trigger firmware (skip events with >20 hits)"
   echo "  --dt-max VALUE      Maximum time window for efficiency (default: -100)"
   echo "  --dt-min VALUE      Minimum time window for efficiency (default: -180)"
   echo "  -h, --help          Display this help message"
   echo ""
   echo "EXAMPLES:"
   echo "  $0 1000 test_run --batch"
   echo "  $0 1000 test_run --batch --dt-max -120 --dt-min -200"
   echo "  $0 1000 test_run --self --dt-max -150 --dt-min -220"
   exit 0
fi

# Check for required arguments
if [ "$#" -lt 2 ]; then
   echo "ERROR: Missing required arguments."
   echo "Usage: $0 <nEvents> <comment> [OPTIONS]"
   echo "Use --help for more information."
   exit 1
fi

# Define variables & paths
scriptDir="$(cd "$(dirname "$0")" && pwd)"
rootDir="$(cd "$scriptDir/.." && pwd)"
rootMacrosDir="$scriptDir/root_macros"

dataDir="$rootDir/data"                # Base data directory
inputDir="$dataDir/input"              # Used for storing raw DCT tmp.strip files that are later used in C++ analysis code (hence why it's called input)
outputDir="$dataDir/output"            # This is where all the ROOT files and plots will be stored, organized by timestamp

nEvents="$1"
comment="$2"

timestamp=$(date +"%Y-%m-%d_%H-%M")
tcl_script="GO.tcl"

# Initialize time window parameters with defaults
dt_max="-100"
dt_min="-180"

# Process flags and optional arguments
vivado_mode=""
firmwareDir="BI_DCT_FW"
use_self_trigger=""

# Parse all arguments starting from the 3rd
for arg in "${@:3}"; do
    if [ "$arg" == "-b" ] || [ "$arg" == "--batch" ]; then
        vivado_mode="-mode batch"
    elif [ "$arg" == "-s" ] || [ "$arg" == "--self" ]; then
        firmwareDir="BI_DCT_FW_selftrigger"
        use_self_trigger="--self"
    elif [ "$arg" == "--dt-max" ]; then
        # Next argument should be the value
        shift
        if [ $# -gt 0 ]; then
            dt_max="$1"
        fi
    elif [ "$arg" == "--dt-min" ]; then
        # Next argument should be the value
        shift
        if [ $# -gt 0 ]; then
            dt_min="$1"
        fi
    fi
done

# Update TCL script with parameters
sed -i "1s|.*|set nEvents ${nEvents}|" "$scriptDir/$tcl_script"
sed -i "2s|.*|set timestr ${timestamp}|" "$scriptDir/$tcl_script"
sed -i "3s|.*|set firmwareDir ${firmwareDir}|" "$scriptDir/$tcl_script"
sed -i "4s|.*|set inputDir ${inputDir}|" "$scriptDir/$tcl_script"

# Create output directories for raw and root files if they don't exist
if [ ! -d "$inputDir" ]; then
   mkdir -p "$inputDir"
fi

if [ ! -d "$outputDir" ]; then
   mkdir -p "$outputDir"
fi

mkdir -p "$inputDir/$timestamp"

# Run Vivado data acquisition script
echo "Start acquisition of $nEvents events..."
vivado $vivado_mode -source "$scriptDir/$tcl_script"

# Process acquired data
cd "$inputDir/$timestamp"

python3 "$scriptDir/merge_ila_files.py"
rm -f tmp_file*

grep -v -e Sample -e Rad iladata.txt | awk -F"," '{print $1,$4,$5}' > tmp.strip

# Compile main C++ analysis executable if not already compiled
cd "$rootDir"
if [ ! -f "$rootDir/bin/analysis" ]; then
    echo "Compiling main analysis executable..."
    make build
fi

# Run analysis executable to process raw hits and generate ROOT file
echo "Processing raw hits and generating ROOT file..."
"$rootDir/bin/analysis" "$inputDir/$timestamp/tmp.strip" "$dt_max" "$dt_min" $use_self_trigger

rootfile_name="${timestamp}_${comment}.root"
cp ./output.root "$outputDir/$rootfile_name"

# Cleanup
rm -f "$rootDir"/vivado*
rm -f "$scriptDir"/vivado*