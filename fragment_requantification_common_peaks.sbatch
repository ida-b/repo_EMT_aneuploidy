#!/bin/bash

#SBATCH -o logs/fragment_requantification_common_peaks_%J.txt
#SBATCH -e logs/fragment_requantification_common_peaks_%J.txt
#SBATCH -J fragment_requantification_common_peaks
#SBATCH --ntasks 1
#SBATCH --mem=128G
#SBATCH -t 06:00:00
#SBATCH -p slim16

# Load conda environment
source /store24/project24/ladcol_023/miniconda3/etc/profile.d/conda.sh
conda activate scanalysis

# Run the specified python script
python fragment_requantification_common_peaks.py