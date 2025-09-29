#!/bin/bash

# Take the directory to be the first positional argument when the bash file is run. 
dir=$1

count=$(find "$dir" -type f | wc -l )

echo "There are $count files in the $dir directory."
