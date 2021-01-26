#!/bin/bash
#SBATCH -p broadwl
#SBATCH --qos=covid-19
#SBATCH --account=covid-19
#SBATCH --ntasks-per-node=28
#SBATCH -N 1
#SBATCH --array=11
#SBATCH --mem=32GB
#SBATCH --exclude=midway2-bigmem05,midway2-0217
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --output=logs/reg_%a_mif_%A.out
#SBATCH --error=logs/reg_%a_mif_%A.err

module load R/4.0.0
Rscript fit_regions.R 28