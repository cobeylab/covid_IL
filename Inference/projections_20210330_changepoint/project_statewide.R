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

deltaT = 0.1

print("Writing output")

compartments = c(`Total infectious`='cases_',
  Inc='cases_new_',
  new_deaths='deaths_',
  ObsDeaths='deaths_det_',
  IH_nonicu='hosp_bed_',
  R='recovered_',
  IC='icu_',
  IH='total_hosp_',
  D='cumulative_total_deaths_',
  S='S',
  ObsICU='icu_obs_',
  ObsHosp='hosp_obs_',
  new_hospitalizations='new_hospitalizations_')

## Output civis format
foreach (r=1:11, .combine='rbind') %do%{
  p = read.csv(sprintf("%s/projections_%s_%s.csv", human_output_dir, model_name, r)) %>%
    mutate(Date=as.Date(date)) %>%
    select(-date) %>%
    mutate(Compartment = as.character(Compartment)) %>%
    filter(Compartment %in% names(compartments)) %>%
    mutate(SimID = as.numeric(factor(SimID)))
  p
} ->alloutputs

alloutputs = alloutputs %>%
    group_by(SimID, Date, Compartment) %>%
    summarize(Cases = sum(Cases)) %>%
    ungroup() %>%
    #mutate(Time=as.numeric(Date - as.Date('2020-01-14'))) %>%
    mutate(Region = 'Illinois')

civis_proj = alloutputs %>% rename(date=Date)

print("Projections processed")

write.csv(civis_proj, sprintf("%s/projections_%s_%s.csv", human_output_dir, model_name, 'illinois'))


foreach (i=1:length(compartments)) %do%{
    comp = names(compartments[i])
    comp_name = compartments[i]
    
    cases = civis_proj %>%
        filter(date >= as.Date('2020-03-13')) %>%
        filter(Compartment == comp) %>%
        group_by(date) %>%
        dplyr::summarize(med=quantile(Cases, 0.5, type=3),
                         lower=quantile(Cases, 0.025, type=3),
                         upper=quantile(Cases, 0.975, type=3)) %>%
        mutate(geography_modeled = 'illinois',
               scenario_name = 'baseline') %>%
        ungroup() %>%
        select(date, geography_modeled, scenario_name, med, lower, upper)
    names(cases) = c('date','geography_modeled','scenario_name',paste0(comp_name, c('median','lower','upper')))
    cases
} -> tempframe
finaldf = tempframe %>% reduce(full_join, by = c("date", "scenario_name", "geography_modeled"))
write.csv(finaldf, sprintf('%s/uchicago_covidregion_%s_%s.csv', human_output_dir, 'illinois', Sys.Date()-1), row.names=F)
