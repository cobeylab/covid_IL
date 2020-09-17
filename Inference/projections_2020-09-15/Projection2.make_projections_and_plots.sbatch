#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 28
#SBATCH -N 1
#SBATCH --mem=256gb

cp ../inference_2020-09-15/final_points.csv ./full.final_points.csv
module load R/3.5.1
Rscript Projection1.calculate_R0.R
Rscript Projection2.make_projections_and_plots.R
Rscript Projection3.plot.R