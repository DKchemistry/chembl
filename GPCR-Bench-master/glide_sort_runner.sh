#!/bin/bash

# Define the base command without the input and output file paths
BASE_COMMAND="$SCHRODINGER/utilities/glide_sort -o"

# Prompt for confirmation function
prompt_confirmation() {
    read -p "Proceed with execution? (y/n): " confirm
    if [[ $confirm =~ ^[yY]$ ]]; then
        return 0 # Proceed
    else
        return 1 # Exit
    fi
}

# Generate and echo commands with dry run
generate_and_echo_commands() {
    find . -type f -name "*_lib.sdf" -exec readlink -f {} \; | while read -r file_path; do
        output_path="${file_path%_lib.sdf}_lib_sorted.sdf"
        cmd="$BASE_COMMAND \"$output_path\" -use_dscore -best_by_title \"$file_path\""
        echo "$cmd"
    done | parallel --dry-run
}

# Execute commands
execute_commands() {
    find . -type f -name "*_lib.sdf" -exec readlink -f {} \; | while read -r file_path; do
        output_path="${file_path%_lib.sdf}_lib_sorted.sdf"
        cmd="$BASE_COMMAND \"$output_path\" -use_dscore -best_by_title \"$file_path\""
        echo "$cmd"
    done | parallel
}

# Main script execution
generate_and_echo_commands
if prompt_confirmation; then
    echo "Executing commands..."
    execute_commands
else
    echo "Execution canceled."
fi

