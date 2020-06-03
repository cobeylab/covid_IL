#!/bin/bash
#SBATCH -p cobey
#SBATCH -n 1
#SBATCH -N 1

module load R/3.5.1
Rscript simulate_and_plot.R ${output_dir} ${measurement_model} ${cases} ${population1} ${population2} ${population3} ${t0} ${tf} ${data_start} ${intervention_start} ${ncores} ${maxjobs} ${ICU_start} ${init_file} ${rproc_file} ${nu_scale_file} ${frac_underreport_file} ${contact_file} ${inference_file} ${t_ref} ${intervention_file} ${population4} ${region_to_test}