library(dplyr)
library(pomp)
library(tidyr)
library(foreach)

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

fracs = fraction_underreported %>% 
    mutate(excess = observed_all_cause - exp_all_cause, 
        excess=ifelse(excess<0,0,excess), 
        unreported_deaths = excess - ll_deaths, 
        unreported_deaths=ifelse(unreported_deaths<0,0,unreported_deaths),
        underreport_frac = unreported_deaths / excess) %>%
    select(time, underreport_frac)
