#!/bin/bash

#SBATCH -o logs/processing_%J.txt
#SBATCH -e logs/processing_%J.txt
#SBATCH -J processing
#SBATCH --ntasks 1
#SBATCH --mem=128G
#SBATCH -t 02:00:00
#SBATCH -p slim16

# Load conda environment
source /store24/project24/ladcol_023/miniconda3/etc/profile.d/conda.sh
conda activate scanalysis

# Run the specified python script
python processing.py