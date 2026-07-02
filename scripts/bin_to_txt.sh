#!/usr/bin/env bash
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Use -h or --help for more information."
    exit 1
}

# Skript to convert .pcapng files into .txt files

# Check for help flag first
if [[ "$@" == *"--help"* ]] || [[ "$@" == *"-h"* ]]; then
    echo "Script usage: $0 [OPTIONS]"
    echo ""
    echo "REQUIRED ARGUMENTS:"
    echo "  -d | --data-dir PATH    Path to the data directory containing 'raw' subdirectory"
    echo ""
    echo "OPTIONS:"
    echo "  -i | --input PATH       Convert only this .pcapng file (can be repeated)"
    echo "  -f | --force            Overwrite the target .txt if it already exists"
    echo "  -h | --help             Display this help message"
    echo ""
    echo "EXAMPLES:"
    echo "  $0 --input data/raw/bin/sample.pcapng"
    echo "  $0 --force"
    echo ""
    exit 0
fi

# Check for required arguments
if [ "$#" -lt 1 ]; then
    echo "ERROR: Missing required arguments."
    usage
fi

# Process CLI arguments
dataDir=""
inputs=()
force=0
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -d|--data-dir)
            dataDir="$2"
            shift 2
            ;;
        -i|--input)
            inputs+=("$2")
            shift 2
            ;;
        -f|--force)
            force=1
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Set paths
rootDir="$(dirname "$(dirname "$(realpath "$0")")")"
binDir="$dataDir/raw/bin"
txtDir="$dataDir/raw/txt"

# Check if the folder exists
if [ ! -d "$binDir" ]; then
  echo "Error: folder not found: $binDir"
  exit 1
fi

# Check if tshark is installed
if ! command -v tshark >/dev/null 2>&1; then
  echo "Error: tshark not found"
  exit 1
fi

# Avoid iterating over a fake file named "*.pcapng" if no files are found
shopt -s nullglob

files=()
if [ ${#inputs[@]} -gt 0 ]; then
    for input in "${inputs[@]}"; do
        if [ -f "$input" ]; then
            files+=("$input")
        elif [ -f "$binDir/$input" ]; then
            files+=("$binDir/$input")
        else
            echo "Missing input file: $input"
        fi
    done
else
    # Collect all .pcapng files in the folder
    files=("$binDir"/*.pcapng)
fi

# Exit if no files are found
if [ ${#files[@]} -eq 0 ]; then
  echo "No .pcapng files found in $binDir"
  exit 0
fi

echo "Converting binary into text files"

# Loop through each file once
for file in "${files[@]}"; do
    # Output text file name
    txtFile="$txtDir/$(basename "${file%.pcapng}.txt")"
    tmpTxtFile="${file%.pcapng}.txt"

    # Check if the txt file already exists and handle the force flag
    if [ $force -eq 1 ] && [ -f "$txtFile" ]; then
        rm "$txtFile"
    fi

    # Remove txt file if already exists
    if [ -f "$tmpTxtFile" ]; then
        rm "$tmpTxtFile"
    fi

    # Check if file still exists before processing
    if [ ! -f "$file" ]; then
        echo "File disappeared before processing: $(basename "$file")"
        continue
    fi

    echo "Converting $(basename "$file")..."

    # Run tshark and extract fields
    if tshark -r "$file" -Y "eth.src == ca:02:03:04:05:06" -T fields -e eth.src -e data.data > "$tmpTxtFile"; then
        echo "Done: $(basename "$tmpTxtFile")"
    else
        echo "Error during conversion of $(basename "$file")"
        rm -f "$tmpTxtFile"
        exit 1
    fi

    # Move the txt file to the txt directory
    mv "$tmpTxtFile" "$txtFile"
done

echo "Conversion completed."