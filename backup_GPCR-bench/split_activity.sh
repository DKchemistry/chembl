#!/bin/bash

echo "Starting to process directories..."

# Loop through directories
for dir in */; do
    echo "Processing directory: $dir"
    
    # Ensure there's a .smi file to process
    smi_files=$(find "$dir" -name "*.smi")
    for smi_file in $smi_files; do
        echo "Found .smi file: $smi_file"

        # Prepare output files
        active_file="${smi_file%.smi}_active.smi"
        inactive_file="${smi_file%.smi}_inactive.smi"
        
        echo "Active file will be: $active_file"
        echo "Inactive file will be: $inactive_file"

        # Initialize or clear files
        > "$active_file"
        > "$inactive_file"

        # Process .smi file
        while IFS= read -r line || [[ -n "$line" ]]; do
            # Skip header or empty lines
            if [[ "$line" == "smiles ID" || -z "$line" ]]; then
                echo "Skipping header or empty line"
                continue
            fi

            # Extract ID from line
            id=$(echo "$line" | awk '{print $NF}') # ID is the last field
            
            # Check if ID exists in active or inactive .id files
            if grep -q "$id" "${dir}"/*actives.id 2>/dev/null; then
                echo "Writing to active: $line"
                echo "$line" >> "$active_file"
            elif grep -q "$id" "${dir}"/*decoys.id 2>/dev/null; then
                echo "Writing to inactive: $line"
                echo "$line" >> "$inactive_file"
            else
                echo "ID $id not found in any .id file"
            fi
        done < "$smi_file"
    done
done

echo "Processing complete."

