#!/bin/bash
usage() {
    echo "Usage: $0 --config <config_file>"
    echo "Use --help for more information."
    exit 1
}

# Script to run the analysis workflow based on a config file

# Check for help flag first
if [[ "$@" == *"--help"* ]] || [[ "$@" == *"-h"* ]]; then
    echo "Script usage: $0 --config <config_file>"
    echo ""
    echo "REQUIRED ARGUMENTS:"
    echo "  --config <config_file>  Path to the YAML config file"
    echo ""
    echo "OPTIONS:"
    echo "  --recompile             Force recompilation of the main analysis executable before running"
    echo "  -h, --help              Display this help message"
    exit 0
fi

# Check for required arguments
if [ "$#" -lt 2 ]; then
    echo "ERROR: Missing required arguments."
    usage
fi

# Process CLI arguments
config_file=""
recompile=false

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --config)
            config_file="$2"
            shift 2
            ;;
        --recompile)
            recompile=true
            shift
            ;;
        --help|-h)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

if [ -z "$config_file" ]; then
    echo "ERROR: --config is required."
    usage
fi

if [ ! -e "$config_file" ]; then
    echo "ERROR: Config '$config_file' does not exist."
    exit 1
fi

rootDir="$(dirname "$(dirname "$(realpath "$0")")")"
echo "Root directory: $rootDir"

# Recompile main C++ analysis executable if not already compiled
cd "$rootDir"
if [ "$recompile" = true ]; then
    echo "Recompilation requested. Cleaning previous builds..."
    make clean
    echo "Recompiling analysis executable..."
    make -j$(nproc)
else
    if [ ! -f "$rootDir/bin/analysis" ]; then
        echo "Analysis executable not found. Compiling..."
        make -j$(nproc)
    else
        echo "Analysis executable already exists. Skipping compilation."
    fi
fi

# Run analysis executable on the config file
"$rootDir/bin/analysis" analyze --config "$config_file"

# Update config file with the location of the summary file for plotting
config_base="$(basename "$config_file")"
config_stem="${config_base%.*}"
summary_root="$rootDir/data/output/${config_stem}_summary.root"

if grep -q '^summary root file:' "$config_file"; then
    sed -i '' "s|^summary root file:.*|summary root file: $summary_root|" "$config_file"
else
    printf '\nsummary root file: %s\n' "$summary_root" >> "$config_file"
fi