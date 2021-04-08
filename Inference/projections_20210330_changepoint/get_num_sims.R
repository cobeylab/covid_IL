## Run mif searches, can be implemented in parallel on the cluster 
## Sylvia Ranjeva, Phil Arevalo
## April 2020
## -----------------------------------------------------

## Load packages and specify files 
library(dplyr)
library(pomp)
library(ggplot2)
library(magrittr)
library(lubridate)
library(MASS)
library(reshape2)
library(tidyr)
library(gridExtra)
library(matrixStats)

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)
# Read in inference functions
source('./input_file_specification_no_age.R')
source(covid_get_path(function_file))
source('./set_up_covariates_and_data_no_age.R')

## Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1
num_cores = as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
reg_test = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


stopifnot(reg_test %in% c(-1, 1:11))

outdf = read.csv(sprintf('agg_mif_%s.csv', model_name))


#####

mle_frame = outdf %>%
  filter(region_to_test == reg_test) %>%
  rename(file = filename) %>%
  rename(loglik = loglik.pf) %>%
  filter(!is.na(loglik)) %>%
  filter(region_to_test == reg_test) %>%
  mutate(probs = exp(loglik - logSumExp(loglik)),
    replicates = rmultinom(1, n_sim, probs)) %>%
  arrange(-loglik)

write.csv(mle_frame, sprintf("%s/scenario_num_sims_%s_%s.csv", output_dir, model_name, reg_test))
