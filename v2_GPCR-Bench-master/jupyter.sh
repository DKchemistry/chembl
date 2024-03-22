#!/bin/bash

# Ask for the title suffix
read -p "Enter the title suffix (e.g., S1PR1): " title_suffix

echo "# title suffix"
echo "title_suffix = \"$title_suffix\""
echo ""

# Initialize variables to ensure they're clear if not set
file_path_sdf_active=""
file_path_sdf_decoy=""
file_path_strain_active=""
file_path_strain_decoy=""

# Function to select files with fzf based on pattern
select_files() {
    local pattern=$1
    local prompt_message=$2
    find . -type f -name "${pattern}" | fzf --prompt="${prompt_message}"
}

# Select .sdf and .csv files separately
echo "# Files we are processing"
file_path_sdf_active=$(select_files "${title_suffix}*active*_lib_sorted.sdf" "Select active SDF file: ")
file_path_sdf_decoy=$(select_files "${title_suffix}*decoy*_lib_sorted.sdf" "Select decoy SDF file: ")

file_path_strain_active=$(select_files "${title_suffix}*active*_lib_sorted.csv" "Select active CSV file: ")
file_path_strain_decoy=$(select_files "${title_suffix}*decoy*_lib_sorted.csv" "Select decoy CSV file: ")

# Output for .sdf files
if [ ! -z "$file_path_sdf_active" ]; then
    echo "file_path_sdf_active=\"$file_path_sdf_active\""
else
    echo "file_path_sdf_active=\"\""
fi

if [ ! -z "$file_path_sdf_decoy" ]; then
    echo "file_path_sdf_decoy=\"$file_path_sdf_decoy\""
else
    echo "file_path_sdf_decoy=\"\""
fi

# Output for .csv files
if [ ! -z "$file_path_strain_active" ]; then
    echo "file_path_strain_active = \"$file_path_strain_active\""
else
    echo "file_path_strain_active = \"\""
fi

if [ ! -z "$file_path_strain_decoy" ]; then
    echo "file_path_strain_decoy = \"$file_path_strain_decoy\""
else
    echo "file_path_strain_decoy = \"\""
fi
