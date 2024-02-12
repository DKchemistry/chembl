#!/bin/bash

# Target directory as the first argument (optional)
target_dir="${1:-.}"

# Navigate to the target directory if specified
cd "$target_dir" || { echo "Failed to change directory to $target_dir"; exit 1; }

# Process files
for file in *; do
  if [[ -f "$file" ]]; then # Ensure it's a file
    # Extract ID using extended regex
    id=$(echo "$file" | grep -o -E '^[A-Z0-9]+')
    if [[ -n "$id" ]]; then
      mkdir -p "$id" && mv "$file" "$id/" || echo "Failed to move $file to $id/"
    else
      echo "No ID found for $file"
    fi
  fi
done

