#!/bin/bash

unique_jobid=1010 # An identifier to keep track of runs
task_id=1 # index of the point to sample, currently goes from 1 to the number of rows in grid_search_params.csv
cases=../grid_search/emr_linelist_fitting_data.csv # File with death and ICU data, DO NOT COMMIT TO REPO
output_dir=results_scratch/  # Directory created if it does not already exist

# Shouldn't need to change anything below this point
nu_scale_file=Data/nu_scaling.csv
frac_underreport_file=Data/frac_underreported.csv
measurement_model=Inference/dmeasure_hosp_non_hosp_ICU.c

# Demographic information
population1=Data/IL_population_north-central.csv
population2=Data/IL_population_central.csv
population3=Data/IL_population_northeast.csv
population4=Data/IL_population_southern.csv

init_file=Inference/initializer_non_hosp_deaths.c
rproc_file=Inference/rprocess_calculate_phis.c
contact_file=Data/IL_4metaregions_symmetric.RData
inference_file=Inference/inference_functions.R
intervention_file=Forecasting/intervention_table.csv

t_ref=2020-01-14
t0=2020-03-01
tf=2020-06-04
data_start=2020-03-16
intervention_start=2020-03-16
ICU_start=2020-04-07

ncores=-1 # only necessary if parallelizing via slurm
maxjobs=-1 # only necessary if parallelizing via slurm

module load R/3.5.1
Rscript pfilter_search_theta.R ${task_id} ${unique_jobid} ${output_dir} ${measurement_model} ${cases} ${population1} ${population2} ${population3} ${t0} ${tf} ${data_start} ${intervention_start} ${ncores} ${maxjobs} ${ICU_start} ${init_file} ${rproc_file} ${nu_scale_file} ${frac_underreport_file} ${contact_file} ${inference_file} ${t_ref} ${intervention_file} ${population4}