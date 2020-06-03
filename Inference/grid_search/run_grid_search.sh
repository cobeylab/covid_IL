#!/bin/bash
partitionspec="-p cobey"

export region_to_test=-1

export output_dir=results_region_${region_to_test}/
export nu_scale_file=Data/nu_scaling.csv
export frac_underreport_file=Data/frac_underreported.csv
export measurement_model=Inference/dmeasure_hosp_non_hosp_ICU.c

# File with death and ICU data
export cases=./emr_linelist_fitting_data.csv

export population1=Data/IL_population_north-central.csv
export population2=Data/IL_population_central.csv
export population3=Data/IL_population_northeast.csv
export population4=Data/IL_population_southern.csv

export init_file=Inference/initializer_non_hosp_deaths.c
export rproc_file=Inference/rprocess_calculate_phis.c
export contact_file=Data/IL_4metaregions_symmetric.RData
export inference_file=Inference/inference_functions.R
export intervention_file=Forecasting/intervention_table.csv

export t_ref=2020-01-14
export t0=2020-03-01
export tf=2020-06-01
export data_start=2020-03-16
export intervention_start=2020-03-16
export ICU_start=2020-04-07

export ncores=28
export maxjobs=400

job1=$(sbatch -p cobey --array=1-16 --export=ALL run_pfilter.sbatch)
job2=$(sbatch -p broadwl --qos=covid-19 --account=covid-19 --array=17-100 --export=ALL run_pfilter.sbatch)
job3=$(sbatch ${partitionspec} --mem=16gb --dependency=${job1##* },${job2##* } --export=ALL run_simulate_and_plot.sbatch)
