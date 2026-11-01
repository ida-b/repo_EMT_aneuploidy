#!/bin/bash

#SBATCH -o logs/test_peak_merging_modes_%J.txt
#SBATCH -e logs/test_peak_merging_modes_%J.txt
#SBATCH -J test_peak_merging_modes_
#SBATCH --ntasks 1
#SBATCH --mem=128G
#SBATCH -t 12:00:00
#SBATCH -p slim16

# Load conda environment
source /store24/project24/ladcol_023/miniconda3/etc/profile.d/conda.sh
conda activate scanalysis

# Run the specified python script
python test_peak_merging_modes.py