#!/bin/bash

output_dir=results_scratch/  # Directory created if it does not already exist

module load R/3.5.1
Rscript aggregate_points.R ${output_dir}