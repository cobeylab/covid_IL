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

## Set working directory
cd /project2/cobey/covid-modeling/rt-pipeline-09-2020/estimates/

## Load R
module load R/4.0.0

## Get latest data
cp ${latest_cli_file} ../../internal_data/cli_admissions_latest.csv 
echo "copied latest data to latest.csv"

## Run the Rt estimation pipeline
Rscript estimate_cli_epinow2.R --var=$SLURM_ARRAY_TASK_ID --midway=TRUE --debug=FALSE --outpath=$OP
echo "Ran epinow2"

