#!/bin/bash

#SBATCH --job-name=[whatever_name_you_want]
#SBATCH --partition=cosmology
#SBATCH --output=[whatever_directory_you_want]
#SBATCH --mail-user=[email]
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH  --cpus-per-task=[number of cpu cores]

# Activate conda environment
source ~/.bashrc
conda activate [your_env]

# Run program.
cd [file_loc]
python -u [.py_file] -n [number of cpu cores]

#sbatch job.sh
