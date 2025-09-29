#!/bin/bash

# Source directory containing simulation result folders (must be the first positional argument when running this bash file)
source_dir=$1

# Destination directory where renamed files will be moved (must be the second positional argument when running this bash file)
destination_dir=$2 

# Loop through each folder in the source directory
for folder in "$source_dir"/*/; do
        folder_name=$(basename "$folder"/)  # Get the name of the current folder
        echo "$folder_name"
        # Check if the specific file exists in the current folder
        if [ -f "$folder/data_log.dat" ] ; then
                # Copy the data_log.dat file to the destination directory with the folder name as the new filename
                cp "$folder/data_log.dat" "$destination_dir/$folder_name.dat"
        fi
done

# Count the number of files in the final directory
count=$(find "$destination_dir" -type f | wc -l )

# Print the number of files that were moved
echo "Moved $count files from $source_dir to $destination_dir"
