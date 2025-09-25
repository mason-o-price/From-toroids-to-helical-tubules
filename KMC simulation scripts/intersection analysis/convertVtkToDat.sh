#!/bin/bash

inputDirectory=$1

for folder in "$inputDirectory"/*/; do
    folderName=$(basename "$folder"/)
    echo "Processing $folderName"

    conversions="$inputDirectory/$folderName/conversions"

    for f in "$conversions"/*.vtk; do
        base=$(basename "$f" .vtk) 
        cp "$f" "$conversions/$base.dat"
    done
done
