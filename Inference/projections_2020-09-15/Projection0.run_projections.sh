#!/bin/bash
partitionspec="-p bigmem2"
covidspec="-p bigmem2"


#job1=$(sbatch ${covidspec} Projection1.run_calculate_R0.sbatch)

job2=$(sbatch --dependency=4470750 Projection2.make_projections_and_plots.sbatch)
#--dependency=${job1##* } 