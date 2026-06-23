#!/bin/bash
usage() {
    echo "Usage: $0 --configs <config_file1> <config_file2> ... --output-dir <output_directory> [--recompile]"
    echo "Use --help for more information."
    exit 1
}

# Check for help flag first
if [[ "$@" == *"--help"* ]] || [[ "$@" == *"-h"* ]]; then
    echo "Script usage: $0 --configs <config_files...> --output-dir <output_directory>"
    echo ""
    echo "REQUIRED ARGUMENTS:"
    echo "  --configs <config_files...>      Space-separated paths to YAML config files"
    echo "  --output-dir <output_directory>  Path to the output directory"
    echo ""
    echo "OPTIONS:"
    echo "  --recompile                      Force recompilation of the main analysis executable before running"
    echo "  -h, --help                       Display this help message"
    exit 0
fi

# Check for minimum required arguments
if [ "$#" -lt 2 ]; then
    echo "ERROR: Missing required arguments."
    usage
fi

# Initialize an array to hold configuration files
recompile=false
config_files=()

# Process CLI arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --configs)
            shift
            # Keep pulling files until hitting another option or running out of arguments
            while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
                config_files+=("$1")
                shift
            done
            ;;
        --output-dir)
            output_directory="$2"
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

# Validate captured config files
if [ ${#config_files[@]} -eq 0 ]; then
    echo "ERROR: At least one configuration file is required via --configs."
    usage
fi

# Validate output directory
if [ -z "$output_directory" ]; then
    echo "ERROR: --output-dir is required."
    usage
fi

# Verify every single passed configuration file exists
for config in "${config_files[@]}"; do
    if [ ! -e "$config" ]; then
        echo "ERROR: Config file '$config' does not exist."
        exit 1
    fi
done

rootDir="$(dirname "$(dirname "$(realpath "$0")")")"
echo "Root directory: $rootDir"

# Recompile main C++ analysis executable if requested, otherwise check if it exists
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

# Construct the exact multi-config argument syntax your C++ binary expects
# e.g., --config path1.yaml --config path2.yaml
cpp_args=()
for config in "${config_files[@]}"; do
    cpp_args+=("--config" "$config")
done

echo "Executing plotter with ${#config_files[@]} configurations..."
# Run analysis executable with the expanded arguments array
"$rootDir/bin/analysis" plotter "${cpp_args[@]}" --output-dir "$output_directory"