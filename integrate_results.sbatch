#!/bin/bash

#SBATCH -o logs/integrate_results_%J.txt
#SBATCH -e logs/integrate_results_%J.txt
#SBATCH -J integrate_results
#SBATCH --ntasks 1
#SBATCH --mem=128G
#SBATCH -t 02:00:00
#SBATCH -p slim16

# Load conda environment
source /store24/project24/ladcol_023/miniconda3/etc/profile.d/conda.sh
conda activate scanalysis

# Run the specified python script
python integrate_results.py