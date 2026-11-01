#!/bin/bash

#SBATCH -o logs/gene_activity_analysis_%J.txt
#SBATCH -e logs/gene_activity_analysis_%J.txt
#SBATCH -J gene_activity_analysis
#SBATCH --ntasks 1
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -p slim16

# Load conda environment
source /store24/project24/ladcol_023/miniconda3/etc/profile.d/conda.sh
conda activate scanalysis

# Run the specified python script
python gene_activity_analysis.py
