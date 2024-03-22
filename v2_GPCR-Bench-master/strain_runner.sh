#!/bin/bash

# Initialize variables
BASE_DIR=""
SCRIPT=""
CONDA_ENV_NAME=""

# Parse named arguments
while getopts ":b:s:e:" opt; do
  case ${opt} in
    b ) BASE_DIR=$OPTARG ;;
    s ) SCRIPT=$OPTARG ;;
    e ) CONDA_ENV_NAME=$OPTARG ;;
    \? ) echo "Invalid option: $OPTARG" 1>&2; exit 1 ;;
    : ) echo "Invalid option: $OPTARG requires an argument" 1>&2; exit 1 ;;
  esac
done
shift $((OPTIND -1))

# Check for required arguments
if [ -z "$BASE_DIR" ] || [ -z "$SCRIPT" ] || [ -z "$CONDA_ENV_NAME" ]; then
    echo "Usage: $0 -b <base-directory> -s <script-path> -e <conda-env-name>"
    exit 1
fi

# Activate the Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV_NAME" || { echo "Failed to activate Conda environment: $CONDA_ENV_NAME"; exit 1; }

# Ensure SCRIPT is an absolute path
if [[ ! "$SCRIPT" = /* ]]; then
    echo "Error: script-path must be an absolute path."
    exit 1
fi

# Generate and execute commands
find "$BASE_DIR" -type f -name "*_lib.sdf" -exec readlink -f {} \; | while read -r file_path; do
    output_path=$(echo "$file_path" | sed "s|/docking/|/strain/|" | sed "s|.sdf$|.csv|")
    cmd="python $SCRIPT -i \"$file_path\" -o \"$output_path\""
    echo "$cmd"
done | parallel --dry-run

# Ask for confirmation to proceed with actual execution
read -p "Proceed with execution? (y/n): " confirm && [[ $confirm == [yY] ]] || exit 1

# If confirmed, execute the commands
echo "Executing commands..."
find "$BASE_DIR" -type f -name "*_lib.sdf" -exec readlink -f {} \; | while read -r file_path; do
    output_path=$(echo "$file_path" | sed "s|/docking/|/strain/|" | sed "s|.sdf$|.csv|")
    cmd="python $SCRIPT -i \"$file_path\" -o \"$output_path\""
    echo "$cmd"
done | parallel

