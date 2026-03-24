#!/bin/bash

# Script to run the full data acquisition and processing workflow for the DCT test setup
if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
   echo "Script usage: $0 <nEvents> <comment> [-batch] [-self]"
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

# Process flags
vivado_mode=""
if [ "$3" == "-batch" ] || [ "$4" == "-batch" ]; then
    vivado_mode="-mode batch"
fi

if [ "$3" == "-self" ] || [ "$4" == "-self" ]; then
    firmwareDir="BI_DCT_FW_selftrigger"
    else
    firmwareDir="BI_DCT_FW"
fi

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
rm tmp_file*

grep -v -e Sample -e Rad iladata.txt | awk -F"," '{print $1,$4,$5}' > tmp.strip

root -l -b -q "$rootMacrosDir/process_raw_hits.cpp"
# root -l -b -q "$rootMacrosDir/basic_plots.C"

rootfile_name="${timestamp}_${comment}.root"
plotfile_name="${timestamp}_${comment}.pdf"
cp ./out.root "$rootfileDir/$rootfile_name"
cp ./out.pdf "$plotfileDir/$plotfile_name"

# Cleanup
cd "$rootDir"
rm "$rootDir"/vivado*