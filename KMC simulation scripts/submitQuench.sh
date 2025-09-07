#!/bin/bash
#SBATCH --account=guest
#SBATCH --partition=guest-compute
#SBATCH --tasks=1
#SBATCH --output=output.txt
#SBATCH --job-name=mason_gcmc
#SBATCH --error=err.out

unset init_file_path

cd $openmesh_src;
$openmesh_bin/Assembly -q -i /home/masonprice/quench_input_422_eb-9_kb1000.0_ensemble37313707.json
