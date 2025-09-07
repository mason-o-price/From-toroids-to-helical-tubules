#!/bin/bash

#SBATCH --account=guest
#SBATCH --partition=guest-compute
#SBATCH --array=1-2450
#SBATCH --error=errQuenchArray.txt
#SBATCH --output=outputQuenchArray.txt
#SBATCH --job-name=gcmcToroid

# declare the location of the input files
json_directory="/home/masonprice/toroids_new_tests/quench_jsons_toroid_mu4_5"

# Find the index of the file in the directory based on the array index for this job
file_indx=$SLURM_ARRAY_TASK_ID

# find the name of the nth json file to run, based on the file index
json_file=$(ls -1 "$json_directory" | sed -n "${file_indx}p")

echo "Processing file; $file_indx; $json_file"

# Unset the variable init_file_path so that we can specify a starting structure outside of the default folder
unset init_file_path

# run the quench simulation for this json file
$openmesh_bin/Assembly -q -i $json_directory/$json_file
