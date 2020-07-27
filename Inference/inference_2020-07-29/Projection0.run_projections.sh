#!/bin/bash
partitionspec="-p cobey"
covidspec="-p broadwl-lc" #"-p broadwl --qos=covid-19 --account=covid-19"


job1=$(sbatch ${covidspec} Projection1.run_calculate_R0.sbatch)

job2=$(sbatch ${covidspec} --dependency=${job1##* } Projection2.make_projections_and_plots.sbatch)