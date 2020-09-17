#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=16gb
#SBATCH --array=1


module load R/3.5.1
Rscript 4.final_aggregate.R ${output_dir}