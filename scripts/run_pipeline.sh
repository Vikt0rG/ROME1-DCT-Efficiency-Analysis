#!/bin/bash
usage() {
    echo "Usage: $0 --config <config_file> --dt-max <dt_max> --dt-min <dt_min> --output-dir <output_directory> [--no-external] [--recompile] [--force]"
    echo "Use --help for more information."
    exit 1
}

# Check for help flag first
if [[ "$@" == *"--help"* ]] || [[ "$@" == *"-h"* ]]; then
    echo "Script usage: $0 --config <config_file> --dt-max <dt_max> --dt-min <dt_min> --output-dir <output_directory>"
    echo ""
    echo "REQUIRED ARGUMENTS:"
    echo "  --config PATH             Path to the YAML config file"
    echo "  --dt-max VALUE            Maximum time difference"
    echo "  --dt-min VALUE            Minimum time difference"
    echo "  --output-dir PATH         Path to the output directory"
    echo ""
    echo "OPTIONS:"
    echo "  --no-external             Disable external trigger usage in analysis"
    echo "  --recompile               Force recompilation of the main analysis executable before running"
    echo "  --force                   Force overwrite of existing output directory"
    echo "  -h, --help                Display this help message"
    exit 0
fi

# Check for required arguments
if [ "$#" -lt 8 ]; then
    echo "ERROR: Missing required arguments."
    usage
fi

# Process CLI arguments
config_file=""
dt_max=""
dt_min=""
output_directory=""
recompile=false
force=""

# Process flags and optional arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --config)
            config_file="$2"
            shift 2
            ;;
        --dt-max)
            dt_max="$2"
            shift 2
            ;;
        --dt-min)
             dt_min="$2"
             shift 2
             ;;
        --output-dir)
            output_directory="$2"
            shift 2
            ;;
        --no-external)
            no_external="--no-external"
            shift
            ;;
        --recompile)
            recompile=true
            shift
            ;;
        --force)
            force="--force"
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

# Validate required arguments
if [ -z "$config_file" ] || [ -z "$dt_max" ] || [ -z "$dt_min" ] || [ -z "$output_directory" ]; then
    echo "ERROR: Missing required arguments."
    usage
fi

# Validate config file existence
if [ ! -f "$config_file" ]; then
    echo "ERROR: Config file '$config_file' does not exist."
    exit 1
fi

# Validate output directory
if [ ! -d "$output_directory" ] ; then
    echo "ERROR: Output directory '$output_directory' does not exist."
    exit 1
fi

# Verify they are integers
if ! [[ "$dt_max" =~ ^-?[0-9]+$ ]] || ! [[ "$dt_min" =~ ^-?[0-9]+$ ]]; then
    echo "ERROR: dt_max and dt_min must be integers."
    exit 1
fi

# Verify max is greater than min
if [ "$dt_max" -le "$dt_min" ]; then
    echo "ERROR: dt_max must be greater than dt_min."
    exit 1
fi

# Get absolute path of the script's directory
rootDir="$(dirname "$(dirname "$(realpath "$0")")")"
echo "Root directory: $rootDir"

# Step 0: Recompile main C++ analysis executable if not already compiled
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

# Step 1: Process data from the config file
echo "Processing data from config file: $config_file"
"$rootDir/scripts/processMeasurements.py" --config "$config_file" --dt-max "$dt_max" --dt-min "$dt_min" --output-dir "$output_directory" $no_external $force

# Step 2: Run the analysis on the processed data
echo "Running analysis on processed data from the config file: $config_file"
"$rootDir/scripts/run_data_analysis.sh" --config "$config_file" --output-dir "$output_directory" $no_external $force

echo "Pipeline completed successfully. Data stored in: $output_directory"