regions="1-11"
partition_spec="--qos=covid-19 --account=covid-19 -p broadwl"
#partition_spec2="-p cobey"
#partition_spec="-p cobey"

export ifr_constraint=0.008
job1=$(sbatch --array=${regions} --export=ALL ${partition_spec} fit_regions.sbatch)
export ifr_constraint=0.009
job2=$(sbatch --array=${regions} --export=ALL ${partition_spec} fit_regions.sbatch)
job3=$(sbatch --dependency=${job1##* },${job2##* } ${partition_spec} aggregate_results.sbatch)
sbatch --dependency=${job3##* } --array=${regions} ${partition_spec} plot_all.sbatch

job6=$(sbatch --dependency=${job3##* } --array=${regions} ${partition_spec} get_num_sims.sbatch)
job7=$(sbatch --dependency=${job6##* } --array=${regions} ${partition_spec} project_regions.sbatch)
job8=$(sbatch --dependency=${job7##* } --array=1 ${partition_spec} project_statewide.sbatch)