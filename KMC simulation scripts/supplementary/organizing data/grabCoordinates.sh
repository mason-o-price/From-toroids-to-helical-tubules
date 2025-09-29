#!/bin/bash

# Source directory containing simulation result folders
source_dir=$1

# Destination directory where renamed files will be moved
destination_dir=$2 

# Create variable to count how many files we have moved
count=0

# Loop through each folder in the source directory
for folder in "$source_dir"/*/; do
        folder_name=$(basename "$folder"/)  # Get the name of the current folder
        # Check if the specific file exists in the current folder
	if [ -f "$folder/final_structure/vertices_final.dat" ] ; then
                # Copy the file to the destination directory with the folder name as the new filename (and convert to a CSV)
                cp "$folder/final_structure/vertices_final.dat" "$destination_dir/${folder_name}_coordinates.csv"
		echo "${folder_name}_coordinates.csv"
		
		# Increment the counting variable
		count=$((count+1))
        fi

done

# Print the number of files that were moved
echo "Moved $count files from $source_dir to $destination_dir"
