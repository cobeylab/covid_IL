#!/bin/bash
#SBATCH -p cobey
#SBATCH --ntasks-per-node=28
#SBATCH -N 1
#SBATCH --array=11
#SBATCH --mem=16GB
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --output=logs/mif_reg_%a_mif_%A.out
#SBATCH --error=logs/mif_reg_%a_mif_%A.err

module load R/4.0.0
Rscript fit_regions.R 56 ${ifr_constraint}