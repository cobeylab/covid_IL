#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem=20gb

module load R/3.5.1
Rscript Projection2.make_projections_and_plots.R