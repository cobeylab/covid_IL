#!/bin/bash
#SBATCH -p broadwl
#SBATCH --qos=covid-19
#SBATCH --account=covid-19
#SBATCH -n 28
#SBATCH -N 1
#SBATCH --array=1-11
#SBATCH --mem=32GB
#SBATCH --exclude=midway2-bigmem05,midway2-0217
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --output=logs/o_%A_%a.out
#SBATCH --error=logs/e_%A_%a.err

module load R/4.0.0
Rscript slice_regions.R ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_JOB_ID} 28 28 ${SLURM_ARRAY_TASK_ID}
