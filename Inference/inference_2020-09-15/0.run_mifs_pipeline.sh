#!/bin/bash
covidspec="-p broadwl --qos=covid-19 --account=covid-19"
partitionspec="-p cobey"

# First, copy over linelist data. Needs to have the emr_deaths column.
#cp /project2/cobey/covid-civis/processed_data/emr_linelist_fitting_data_restore.csv ./emr_linelist_fitting_data.csv

# Then, run script to process public linelist data
#module load Anaconda3
#python util.process_county_data.py


export ncores=28

export output_dir='results_firstmif/'

export maxjobs=28
export region_to_test=1
job7=$(sbatch ${partitionspec} --array=1-4 --export=ALL 3.run_mif.sbatch)
export region_to_test=2
job8=$(sbatch ${partitionspec} --array=1-4 --export=ALL 3.run_mif.sbatch)
export region_to_test=3
job9=$(sbatch ${partitionspec} --array=1-4 --export=ALL 3.run_mif.sbatch)
export region_to_test=4
job10=$(sbatch ${partitionspec} --array=1-4 --export=ALL 3.run_mif.sbatch)
export region_to_test=5
job11=$(sbatch ${covidspec} --array=1-4 --export=ALL 3.run_mif.sbatch)
export region_to_test=6
job12=$(sbatch ${covidspec} --array=1-4 --export=ALL 3.run_mif.sbatch)


# Run final aggregation
job13=$(sbatch ${partitionspec} --dependency=${job7##* },${job8##* },${job9##* },${job10##* },${job11##* },${job12##* } --export ALL 4.run_final_aggregate.sbatch)

# Plots MLE
sbatch ${partitionspec} --dependency=${job13##* } --export ALL 5.run_plot_mle.sbatch


# Compare MLE

sbatch ${partitionspec} --dependency=${job13##* } --export ALL util.run_compare_mle.sbatch
#--dependency=${job7##* },${job8##* },${job9##* },${job10##* },${job11##* }