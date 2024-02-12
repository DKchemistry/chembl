#!/bin/bash

# Function to display help/documentation
show_help() {
    echo "Usage: ligprep_runner.sh [OPTION]"
    echo "Navigate to directories containing ligprep .sh scripts and optionally execute them."
    echo ""
    echo "Options:"
    echo "  -h          Display this help and exit"
    echo "  -r          Actually execute the scripts; without this, the script runs in dry mode"
    echo ""
    echo "In the default dry mode, the script will navigate to directories and echo the present working directory, without performing any actions."
}

DRY_RUN=1

# Parse options
while getopts "hr" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        r )
            DRY_RUN=0
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            show_help
            exit 1
            ;;
    esac
done

sh_files=$(find . -type f -path "*/stereochemistry/*_ligprep.sh")

if [ -z "$sh_files" ]; then
    echo "No .sh files found to process."
    exit 0
fi

# Process .sh files
echo "$sh_files" | while read sh_file; do
    dir=$(dirname "$sh_file")
    echo "Navigating to $dir"
    cd "$dir" || exit
    echo "Current directory: $(pwd)"
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "Dry run: Would execute $(basename "$sh_file")"
    else
        echo "Executing $(basename "$sh_file")"
        bash "$(basename "$sh_file")"
    fi
    cd - > /dev/null
done
