#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=24GB
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --output=logs/o_%A_%a.out
#SBATCH --error=logs/e_%A_%a.err

module load R/3.5.1
Rscript extract_median_hospital_admissions.R
