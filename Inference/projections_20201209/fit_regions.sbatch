#!/bin/bash
#SBATCH -p broadwl
#SBATCH --account=covid-19
#SBATCH --qos=covid-19
#SBATCH -n 28
#SBATCH -N 1
#SBATCH --array=1
#SBATCH --mem=32GB
#SBATCH --exclude=midway2-bigmem05,midway2-0217
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --output=logs/o_%A_%a.out
#SBATCH --error=logs/e_%A_%a.err

module load R/3.5.1
Rscript fit_regions.R ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_JOB_ID} 28 100 ${region_to_test}
