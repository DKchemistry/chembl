#!/bin/bash

# This script counts occurrences of '$$$$' in input and output .sdf files
# and prints the counts for easy comparison.

find . -type f -name "*_lib.sdf" | while read file; do
  input_file="$file"
  output_file="${file%_lib.sdf}_lib_sorted.sdf"
  
  # Check if the output file exists before trying to count
  if [ -f "$output_file" ]; then
    echo "Input ($input_file): $(grep -c '$$$$' "$input_file")"
    echo "Output ($output_file): $(grep -c '$$$$' "$output_file")"
  else
    echo "Output file $output_file does not exist."
  fi
done

