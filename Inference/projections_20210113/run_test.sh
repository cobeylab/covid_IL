
export region_to_test=11
job11=$(sbatch --array=1 --export=ALL fit_regions.sbatch)
job12=$(sbatch --dependency=${job11##* } aggregate.sbatch)
sbatch --dependency=${job12##* } --array=11 plot_all.sbatch
#job13=$(sbatch --dependency=${job12##* } --array=11 slice_regions.sbatch)
#job14=$(sbatch --dependency=${job13##* } --array=11 project_all.sbatch)