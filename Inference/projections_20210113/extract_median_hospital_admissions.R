library(dplyr)
library(pomp)
library(ggplot2)
library(magrittr)
library(lubridate)
library(MASS)
library(reshape2)
library(tidyr)
library(foreach)
library(gridExtra)
library(purrr)

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


# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

print("Writing output")

foreach (r=1:11, .combine='rbind') %do%{
  p = read.csv(sprintf("%s/projections_%s_%s.csv", human_output_dir, model_name, r))
  #p = read.csv(sprintf("take2_projection_outputs_beta_hfr_mu_gamma_constphi_cdc_outputs/projections_beta_hfr_mu_gamma_constphi_cdc_%s.csv",  r))
  p
} ->alloutputs

alloutputs = alloutputs %>%
    mutate(Date=as.Date(date)) %>%
    select(-date) %>%
    filter(Compartment == 'new_hospitalizations') %>%
    group_by(Date, Compartment, Region) %>%
    summarize(Cases = quantile(Cases, 0.5, type=3)) %>%
    ungroup() %>%
    rename(covid_region = Region, cli=Cases, date=Date) %>%
    select(-Compartment)

write.csv(alloutputs, sprintf('cli_from_model_%s.csv', Sys.Date()), row.names=F)
