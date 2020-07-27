#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 28
#SBATCH -N 1
#SBATCH --mem=16gb

module load R/3.5.1
Rscript Projection1.calculate_R0.R