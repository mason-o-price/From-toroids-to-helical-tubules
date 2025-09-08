#!/bin/bash
#SBATCH --account=guest
#SBATCH --partition=guest-compute
#SBATCH --tasks=1
#SBATCH --output=output.txt
#SBATCH --job-name=mason_gcmc
#SBATCH --error=err.out

json_file = ""

echo "$json_file"
cd $openmesh_src;
$openmesh_bin/Assembly -e -i $json_file

