## Set outpath
export OP="./epinow2_cli_estimates/2020-11-24_cli_v1.2.1/"
export latest_cli_file="../../internal_data/cli_admissions_2020-11-24.csv"
job1=$(sbatch --export=ALL Midway_run_epinow2_cli.sbatch)
sbatch --dependency=${job1##* } --export=ALL Midway_consolidate.sbatch