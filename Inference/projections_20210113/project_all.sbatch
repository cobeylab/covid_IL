#!/bin/bash
#SBATCH -p broadwl
#SBATCH --account=covid-19
#SBATCH --qos=covid-19
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=1-11
#SBATCH --mem=24GB
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --output=logs/o_%A_%a.out
#SBATCH --error=logs/e_%A_%a.err

module load R/4.0.0
Rscript project_regions.R ${SLURM_ARRAY_TASK_ID}
