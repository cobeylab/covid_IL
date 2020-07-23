#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 28
#SBATCH -N 1
#SBATCH --array=1

module load R/3.5.1
Rscript 3.pfilter_around_mle.R ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_JOB_ID} ${ncores} ${maxjobs} ${output_dir}