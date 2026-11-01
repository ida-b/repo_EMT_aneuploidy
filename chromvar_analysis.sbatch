#!/bin/bash

#SBATCH -o logs/chromvar_analysis_%J.txt
#SBATCH -e logs/chromvar_analysis_%J.txt
#SBATCH -J chromvar_analysis
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem=128G
#SBATCH -t 08:00:00
#SBATCH -p slim16

source /store24/project24/ladcol_023/miniconda3/etc/profile.d/conda.sh

# Step 1: export ATAC counts to .mtx/.bed/.csv (Python env)
conda activate scanalysis
python export_for_chromvar.py

# Step 2: chromVAR TF-motif deviations + differential tests (R env)
conda activate caf_atac_r
Rscript chromvar_analysis.R
