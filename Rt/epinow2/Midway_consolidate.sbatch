#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --array=1-12
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH --qos=covid-19
#SBATCH --account=covid-19 
#SBATCH -p broadwl
#SBATCH --job-name=v1.2.1_exact

## Load R
module load R/4.0.0
Rscript summarize_epinow2.R --outpath=$OP

