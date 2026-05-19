#!/usr/bin/env bash
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Use --help for more information."
    exit 1
}

# Skript to convert .pcapng files into .txt files

# Check for help flag first
if [[ "$@" == *"--help"* ]]; then
    echo "Script usage: $0 [OPTIONS]"
    echo ""
    echo "OPTIONS:"
    echo "  --input FILE        Convert only this .pcapng file (can be repeated)"
    echo "  --force             Overwrite the target .txt if it already exists"
    echo ""
    echo "EXAMPLES:"
    echo "  $0 --input data/raw/bin/sample.pcapng"
    echo "  $0 --force"
    echo ""
    exit 0
fi

rootDir="$(dirname "$(dirname "$(realpath "$0")")")"

# Path to the folder containing .pcapng files
binDir="$rootDir/data/raw/bin"
txtDir="$rootDir/data/raw/txt"

inputs=()
force=0
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --input)
            inputs+=("$2")
            shift 2
            ;;
        --force)
            force=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

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