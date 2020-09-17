#!/bin/bash
covidspec="-p broadwl-lc"
partitionspec="-p broadwl --qos=covid-19 --account=covid-19"

# First, copy over linelist data. Needs to have the emr_deaths column.
#cp /project2/cobey/covid-civis/processed_data/emr_linelist_fitting_data.csv ./

# Then, run script to process public linelist data
#module load Anaconda3
#python util.process_county_data.py


export ncores=28

export output_dir='region_3_test/'

export maxjobs=28
export region_to_test=3
job9=$(sbatch ${covidspec} --array=1-4 --export=ALL 3.run_mif.sbatch)

# Run final aggregation
job12=$(sbatch ${partitionspec} --dependency=${job9##* } --export ALL 4.run_final_aggregate.sbatch)

# Plots MLE
sbatch ${covidspec} --export ALL --dependency=${job12##* } 5.run_plot_mle.sbatch

