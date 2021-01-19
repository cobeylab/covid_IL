#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 28
#SBATCH -N 1
#SBATCH --array=1-11
#SBATCH --mem=200GB
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --output=logs/o_%A_%a.out
#SBATCH --error=logs/e_%A_%a.err

module load R/3.5.1
Rscript simulate_with_pfilter.R ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_JOB_ID} 15 200 ${SLURM_ARRAY_TASK_ID}
