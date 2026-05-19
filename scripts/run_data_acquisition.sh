#!/bin/bash
usage() {
    echo "Usage: $0 <nEvents> <comment> [OPTIONS]"
    echo "Use --help for more information."
    exit 1
}

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
    usage
fi

# Define variables & paths
scriptDir="$(cd "$(dirname "$0")" && pwd)"
rootDir="$(cd "$scriptDir/.." && pwd)"
rootMacrosDir="$scriptDir/root_macros"

dataDir="$rootDir/data"                # Base data directory
rawDir="$dataDir/raw"                  # Used for storing raw DCT <name>.txt files that are later used in C++ analysis code
outputDir="$dataDir/output"            # This is where all the ROOT files and plots will be stored, organized by timestamp

nEvents="$1"
comment="$2"

timestamp=$(date +"%Y-%m-%d_%H-%M")
outputFileName="${timestamp}_${comment}"
tcl_script="GO.tcl"

# Initialize time window parameters with defaults
dt_max="-100"
dt_min="-180"

# Process flags and optional arguments
vivado_mode=""
firmwareDir="BI_DCT_FW"
use_self_trigger=""

# Parse all arguments starting from the 3rd
shift 2
while [ $# -gt 0 ]; do
    case "$1" in
        -b|--batch)
            vivado_mode="-mode batch"
            shift
            ;;
        -s|--self)
            firmwareDir="BI_DCT_FW_selftrigger"
            use_self_trigger="--self"
            shift
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Switch into the root project directory to ensure relative paths work correctly
cd "$rootDir"

# Update TCL script with parameters
sed -i "1s|.*|set nEvents ${nEvents}|" "$scriptDir/$tcl_script"
sed -i "2s|.*|set timestr ${timestamp}|" "$scriptDir/$tcl_script"
sed -i "3s|.*|set firmwareDir ${firmwareDir}|" "$scriptDir/$tcl_script"
sed -i "4s|.*|set rawDir ${rawDir}|" "$scriptDir/$tcl_script"

# Create output directories for raw and root files if they don't exist
if [ ! -d "$rawDir" ]; then
   mkdir -p "$rawDir"
fi

if [ ! -d "$outputDir" ]; then
   mkdir -p "$outputDir"
fi

# mkdir -p "$rawDir/$timestamp"
# cd "$rawDir/$timestamp"
cd "$rawDir"

# Run Vivado data acquisition script
echo "Start acquisition of $nEvents events..."
vivado $vivado_mode -source "$scriptDir/$tcl_script"

# Merge files, skip headers, and filter in one pass
( for f in tmp_file*; do tail -n +2 "$f"; done ) | \
  grep -v -e Sample -e Rad | \
  awk -F"," '{print $1,$4,$5}' > "$outputFileName.txt"

rm -f tmp_file*

# Cleanup
rm -f "$rootDir"/vivado*
rm -f "$scriptDir"/vivado*