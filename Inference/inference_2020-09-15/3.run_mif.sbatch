#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 28
#SBATCH -N 1
#SBATCH --array=1

module load R/3.5.1
Rscript 3.mif_single.R ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_JOB_ID} ${ncores} ${maxjobs} ${region_to_test} ${output_dir}