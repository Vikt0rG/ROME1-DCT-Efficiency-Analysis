#!/bin/bash
usage() {
    echo "Usage: $0 -c|--config <config_file> --dt-max <dt_max> --dt-min <dt_min> -d|--data-dir <data_directory> [--no-external] [-r|--recompile] [-f|--force]"
    echo "Use -h or --help for more information."
    exit 1
}

# Check for help flag first
if [[ "$@" == *"--help"* ]] || [[ "$@" == *"-h"* ]]; then
    echo "Script usage: $0 -c | --config <config_file> --dt-max <dt_max> --dt-min <dt_min> -d | --data-dir <data_directory>"
    echo ""
    echo "REQUIRED ARGUMENTS:"
    echo "  -c | --config PATH        Path to the YAML config file"
    echo "  --dt-max VALUE            Maximum time difference"
    echo "  --dt-min VALUE            Minimum time difference"
    echo "  -d | --data-dir PATH      Path to the data (where raw/txt/root subdirectories are located) directory"
    echo ""
    echo "OPTIONS:"
    echo "  -p | --plotter            Run the data plotter after analysis"
    echo "  --no-external             Disable external trigger usage in analysis"
    echo "  -r | --recompile          Force recompilation of the main analysis executable before running"
    echo "  -f | --force              Force overwrite of existing output directory"
    echo "  -h | --help               Display this help message"
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
data_directory=""
plotter=false
recompile=false
force=""
# Process flags and optional arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -c|--config)
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
        -d|--data-dir)
            data_directory="$2"
            shift 2
            ;;
        -p | --plotter)
            plotter=true
            shift
            ;;
        --no-external)
            no_external="--no-external"
            shift
            ;;
        -r|--recompile)
            recompile=true
            shift
            ;;
        -f|--force)
            force="--force"
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

# Validate required arguments
if [ -z "$config_file" ] || [ -z "$dt_max" ] || [ -z "$dt_min" ] || [ -z "$data_directory" ]; then
    echo "ERROR: Missing one of the required arguments."
    usage
fi

# Validate config file existence
if [ ! -f "$config_file" ]; then
    echo "ERROR: Config file '$config_file' does not exist."
    exit 1
fi

# Validate data directory
if [ ! -d "$data_directory" ] ; then
    echo "ERROR: Data directory '$data_directory' does not exist."
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

# Remove trailing slash from data directory if present
data_directory="${data_directory%/}"

# Setup and validate paths
rootDir="$(dirname "$(dirname "$(realpath "$0")")")"
output_directory="$data_directory/output"
echo ""
echo "Root directory: $rootDir"
echo "Output directory: $output_directory"

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

# Set and export variables for the pipeline so all the child scripts can access them
export IN_PIPELINE=true
export ROOT_DIR="$rootDir"

# Step 1: Process data from the config file (convert binary to txt and process txt to produce processed root files)
echo ""
echo "Processing data from config file: $config_file"
python3 "$rootDir/scripts/processMeasurements.py" "--config" "$config_file" "--dt-max" "$dt_max" "--dt-min" "$dt_min" "--data-dir" "$data_directory" $no_external $force

# Step 2: Run the analysis on the processed data
echo ""
echo "Running analysis on processed data from the config file: $config_file"
bash "$rootDir/scripts/run_data_analysis.sh" --config "$config_file" --output-dir "$output_directory"

# Step 3 [Optional]: Run plotter on produced summary file for this config file and remove the summary file afterwards
echo ""
if [ "$plotter" = true ]; then
    echo "Running plotter on the produced summary file for the config file: $config_file"
    bash "$rootDir/scripts/run_data_plotter.sh" --config "$config_file" --output-dir "$output_directory"
fi
echo "Pipeline completed successfully. Data stored in: $data_directory"