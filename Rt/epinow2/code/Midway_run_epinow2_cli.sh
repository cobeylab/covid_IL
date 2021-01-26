#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --array=1-12
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=parevalo
#SBATCH -p cobey
#SBATCH --job-name=v1.2.1_exact

## #SBATCH --output=midway/o_%A_%a.out
## #SBATCH	--error=midway/e_%A_%a.err
## SBATCH --qos=covid-19
## SBATCH --account=covid-19 -p broadwl
## SBATCH --partition=broadwl
## Set outpath
OP="../epinow2_cli_estimates/2021-01-26_test/"
latest_cli_file="/project2/cobey/covid-modeling/rt-pipeline-09-2020/data/cli_admissions_2021-01-21.csv"


## Set working directory
cd /project2/cobey/covid-modeling/rt-pipeline-09-2020/estimates/


## Load R
module load R/4.0.0

## Get latest data
##  !!!!! YOU MUST UPDATE THE DATE IN THIS LINE WITH EACH RUN !!!!!!!!!
cp ${latest_cli_file} ../data/cli_admissions_latest.csv 
echo "copied latest data to latest.csv"

## Run the Rt estimation pipeline
Rscript estimate_cli_epinow2.R --var=$SLURM_ARRAY_TASK_ID --midway=TRUE --debug=FALSE --outpath=$OP
echo "Ran epinow2"

## Reformat the outputs, and copy them into a single .csv
Rscript summarize_epinow2.R

