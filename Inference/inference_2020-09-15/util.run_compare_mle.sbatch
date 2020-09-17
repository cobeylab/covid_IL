#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 6
#SBATCH --mem=32gb
#SBATCH -N 1

module load R/3.5.1
Rscript util.compare_mles.R