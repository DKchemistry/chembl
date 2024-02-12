#!/bin/bash

# If the first argument is -h or --help, display the help message and exit
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo "Usage: $0 {directory} [-run]"
    echo
    echo "This script is used to enumerate the stereoisomers of active and decoy molecules in the subdirectory of the given directory."
    echo "The script expects the directory to process as an argument. You can pass . as the argument to process the current directory."
    echo "It uses the RDKitEnumerateStereoisomers.py script to perform the enumeration with '-d no' and '-m UnassignedOnly' options."
    echo "The script will create a 'stereochemistry' folder within each subdirectory if it does not exist."
    echo "The detection of "active" and "decoy" files is based on the presence of '_active.smi' and '_decoy.smi' files in the subdirectory."
    echo "If the -run flag is provided, the script will execute the commands. Otherwise, it will only print them (debugging)."
    echo
    exit 0
fi

# Get the list of subdirectories
subdirs=$(find "$1" -mindepth 1 -maxdepth 1 -type d)

# Iterate over each subdirectory
for subdir in $subdirs; do
    # Check if a "stereochemistry" folder exists within the subdirectory
    if [ ! -d "$subdir/stereochemistry" ]; then
        mkdir "$subdir/stereochemistry"
    fi

    # Get the list of active and decoy files
    active_files=$(find "$subdir" -maxdepth 1 -type f -name "*_active.smi")
    decoy_files=$(find "$subdir" -maxdepth 1 -type f -name "*_decoy.smi")

    # Iterate over each pair of active and decoy files
    for active_file in $active_files; do
        decoy_file=$(echo "$active_file" | sed 's/_active/_decoy/')
        output_file_active=$(basename "$active_file" .smi)_sc.smi
        output_file_decoy=$(basename "$decoy_file" .smi)_sc.smi

        # Determine the working directory for this subdirectory
        working_dir=$(realpath "$subdir/stereochemistry")

        # Get the absolute paths of input files
        active_file=$(realpath "$active_file")
        decoy_file=$(realpath "$decoy_file")

        # Define the output paths based on the realpath of subdirectory
        output_file_active="$working_dir/$(basename "$active_file" .smi)_sc.smi"
        output_file_decoy="$working_dir/$(basename "$decoy_file" .smi)_sc.smi"

        # Run the Python script on the input files and specify the output files
        command="$MCT/RDKitEnumerateStereoisomers.py -d no -m UnassignedOnly --overwrite \
            --workingdir \"$working_dir\" -i \"$active_file\" -o \"$output_file_active\""

        # Add logging command
        command+=" > \"$working_dir/$(basename "$active_file" .smi)_sc.log\" 2>&1"

        # Execute the command if -run flag is provided
        if [ "$2" = "-run" ]; then
            eval $command
        else
            echo $command
        fi

        command="$MCT/RDKitEnumerateStereoisomers.py -d no -m UnassignedOnly --overwrite \
            --workingdir \"$working_dir\" -i \"$decoy_file\" -o \"$output_file_decoy\""

        # Add logging command
        command+=" > \"$working_dir/$(basename "$decoy_file" .smi)_sc.log\" 2>&1"

        # Execute the command if -run flag is provided
        if [ "$2" = "-run" ]; then
            eval $command
        else
            echo $command
        fi
    done
done
