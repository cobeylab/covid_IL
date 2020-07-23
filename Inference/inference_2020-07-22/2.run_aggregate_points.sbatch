#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 28
#SBATCH -N 1
#SBATCH --array=1

module load R/3.5.1
Rscript 2.aggregate_points.R ${output_dir}