#!/bin/bash

# Source directory containing simulation result folders
source_dir=$1

# Create variable to count how many files we have moved
count=0

# Loop through each folder in the source directory
for folder in "$source_dir"/*/; do
        folder_name=$(basename "$folder"/)  # Get the name of the current folder
        echo "Processing $folder_name"

	cd $openmesh_src
	$openmesh_bin/Assembly -c "$source_dir/$folder_name/snapshots.om" -f "$source_dir/$folder_name/conversions/"

	count=$((count+1))
done

# Print the number of files that were moved
echo "Processed $count simulations from $source_dir"
