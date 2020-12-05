export region_to_test=1
job1=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=2
job2=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=3
job3=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=4
job4=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=5
job5=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=6
job6=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=7
job7=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=8
job8=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=9
job9=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=10
job10=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

export region_to_test=11
job11=$(sbatch --array=1 --export=ALL fit_regions.sbatch)

job12=$(sbatch --dependency=${job1##* },${job2##* },${job3##* },${job4##* },${job5##* },${job6##* },${job7##* },${job8##* },${job9##* },${job10##* },${job11##* } aggregate.sbatch)

sbatch --dependency=${job12##* } plot_all.sbatch

job13=$(sbatch --dependency=${job12##* } slice_regions.sbatch)

job14=$(sbatch --dependency=${job13##* } project_all.sbatch)

sbatch --dependency=${job14##* } project_state.sbatch