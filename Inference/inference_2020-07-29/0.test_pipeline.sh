#!/bin/bash
partitionspec="-p cobey"
covidspec="-p broadwl-lc" #"-p broadwl --qos=covid-19 --account=covid-19"

export ncores=4
export maxjobs=4
export output_dir='mifs_testing/'

export region_to_test=1
job1=$(sbatch ${covidspec} --array=1 --export=ALL 1.run_mif.sbatch)

export region_to_test=2
job2=$(sbatch ${covidspec} --array=1 --export=ALL 1.run_mif.sbatch)

export region_to_test=3
job3=$(sbatch ${covidspec} --array=1 --export=ALL 1.run_mif.sbatch)

export region_to_test=4
job4=$(sbatch ${covidspec} --array=1 --export=ALL 1.run_mif.sbatch)

export region_to_test=5
job5=$(sbatch ${covidspec} --array=1 --export=ALL 1.run_mif.sbatch)

job6=$(sbatch ${covidspec} --dependency=${job1##* },${job2##* },${job3##* },${job4##* },${job5##* } --export=ALL 2.run_aggregate_points.sbatch)

export ncores=4
export maxjobs=4
job7=$(sbatch ${covidspec} --dependency=${job6##* } --array=1 --export=ALL 3.run_pfilter_around_mle.sbatch)

job10=$(sbatch ${covidspec}  --dependency=${job7##* } --export ALL 4.run_final_aggregate.sbatch)

sbatch ${covidspec} --export ALL --dependency=${job10##* } 5.run_plot_mle.sbatch
#