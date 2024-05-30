#!/bin/bash

# Define the base directory
base_dir="/Users/lkv206/work/to_do_projects/chembl_ligands/v3_GPCR-Bench-master"

# Navigate to the base directory
cd "$base_dir" || exit

# Loop through each directory and create the freeform_strain directory
for dir in [A-Z0-9]*/; do
    mkdir -p "$dir/freeform_strain"
done

echo "Directories created successfully."
