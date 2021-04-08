#!/bin/bash
#SBATCH -p cobey
#SBATCH --ntasks-per-node 8
#SBATCH -N 1
#SBATCH --array=1-11
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --output=logs/nsim_o_%A_%a.out
#SBATCH --error=logs/nsim_e_%A_%a.err

module load R/4.0.0
Rscript get_num_sims.R