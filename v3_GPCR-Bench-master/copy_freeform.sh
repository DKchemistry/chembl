#!/bin/bash

# Define the source and destination base directories
source_dir="/Users/lkv206/work/to_do_projects/chembl_ligands/stasa/strain_b"
destination_base_dir="/Users/lkv206/work/to_do_projects/chembl_ligands/v3_GPCR-Bench-master"

# Loop through each file in the source directory
for file in "$source_dir"/*_f.csv; do
    # Extract the directory name from the file name (e.g., 5HT1B from 5HT1B_active_f.csv)
    dir_name=$(basename "$file" | cut -d'_' -f1)
    
    # Define the destination directory
    dest="$destination_base_dir/$dir_name/freeform_strain/"
    
    # Check if the destination directory exists
    if [ -d "$dest" ]; then
        # Copy the file to the destination directory, retaining the original file name
        cp "$file" "$dest"
        echo "Copied $(basename "$file") to $dest"
    else
        echo "Directory $dest does not exist."
    fi
done
