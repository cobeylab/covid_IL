## Aggregate epinow2 estimates for each covid region, copy into one data frame and then save to the ../figs/... directory
library(readr)
library(rlang)
library(dplyr)

## Source function used to summarise outputs
source('../code/summarise-epinow2-estimates.R')

## Read in job information saved from script used to generate estimates
runpar_file = sprintf('run_params_cli_%s.rds', as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID")))
runpars <- read_rds(runpar_file)

## Summarise the estimates
summarise_all_estimates(path = runpars$outpath, 
    dt = as.Date(runpars$dt), 
    tooday = Sys.Date())
