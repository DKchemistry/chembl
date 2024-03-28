#!/bin/bash

# Find directories named 'docking' and store the result in a text file
find . -name "docking" -type d > docking_dirs.txt

# Read the text file line by line
while IFS= read -r dir; do
  # Find .sh files in the current 'docking' directory
  find "$dir" -name "*.sh" > sh_files.txt
  
  # Read the list of .sh files and execute them
  while IFS= read -r sh_file; do
    # Extract directory path
    dir_path=$(dirname "$sh_file")
    
    # Change to the script's directory
    cd "$dir_path" || exit
    
    # Make sure the script is executable
    chmod +x "$(basename "$sh_file")"
    
    # Execute the script
    echo "Executing $sh_file"
    "./$(basename "$sh_file")"
    
    # Change back to the original directory
    cd - || exit
  done < sh_files.txt
done < docking_dirs.txt

# Cleanup
rm docking_dirs.txt sh_files.txt

