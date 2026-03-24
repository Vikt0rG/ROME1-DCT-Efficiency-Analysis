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
rootMacrosDir="$rootDir/root_macros"

outputDir="$rootDir/data/output"
rawOutputDir="$outputDir/DCT_raw_files"
rootfileDir="$outputDir/DCT_root_files"
plotfileDir="$outputDir/DCT_plots"

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
sed -i "1s/.*/set nEvents ${nEvents}/" "$scriptDir/$tcl_script"
sed -i "2s/.*/set timestr ${timestamp}/" "$scriptDir/$tcl_script"
sed -i "3s/.*/set firmwareDir ${firmwareDir}/" "$scriptDir/$tcl_script"
sed -i "4s/.*/set outputDir ${outputDir}/" "$scriptDir/$tcl_script"

# Create output directories if they don't exist
if [ ! -d "$outputDir" ]; then
   mkdir -p "$outputDir"
fi

if [ ! -d "$rawOutputDir" ]; then
   mkdir -p "$rawOutputDir"
fi

if [ ! -d "$rootfileDir" ]; then
   mkdir -p "$rootfileDir"
fi

if [ ! -d "$plotfileDir" ]; then
   mkdir -p "$plotfileDir"
fi

mkdir -p "$outputDir/$timestamp"

# Run Vivado data acquisition script
echo "Start acquisition of $nEvents events..."
vivado $vivado_mode -source "$scriptDir/$tcl_script"

# Process acquired data
cd "$outputDir/$timestamp"

python3 "$scriptDir/merge_ila_files.py"
rm -f tmp_file*

grep -v -e Sample -e Rad iladata.txt | awk -F"," '{print $1,$4,$5}' > tmp.strip

# Compile process_raw_hits executable if not already compiled
if [ ! -f "$rootDir/bin/process_raw_hits" ]; then
    echo "Compiling process_raw_hits executable..."
    make -C "$rootDir" process_raw_hits
fi

# Run the process_raw_hits executable with parameters
echo "Running process_raw_hits with dt_max=$dt_max, dt_min=$dt_min, and self_trigger=$use_self_trigger..."
"$rootDir/bin/process_raw_hits" --dt-max "$dt_max" --dt-min "$dt_min" $use_self_trigger

rootfile_name="${timestamp}_${comment}.root"
plotfile_name="${timestamp}_${comment}.pdf"
cp ./out.root "$rootfileDir/$rootfile_name"
cp ./out.pdf "$plotfileDir/$plotfile_name"

# Cleanup
cd "$rootDir"
rm -f "$rootDir"/vivado*