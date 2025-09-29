#!/bin/bash

#SBATCH --account=guest
#SBATCH --partition=guest-compute
#SBATCH --array=1-2450
#SBATCH --error=errArray.txt
#SBATCH --output=outputArray.txt
#SBATCH --job-name=gcmcToroid

# declare the location of the input files
json_directory=""

# Find the index of the file in the directory based on the array index for this job
file_indx=$SLURM_ARRAY_TASK_ID

# find the name of the nth json file to run, based on the file index
json_file=$(ls -1 "$json_directory" | sed -n "${file_indx}p")

echo "Processing file; $file_indx; $json_file"

# run the simulation for this json file
$openmesh_bin/Assembly -e -i $json_directory/$json_file
