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

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

source('./input_file_specification_11regions.R')
default_par_file = './final_mle_pars.csv'

# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1

output_dir='./'
region_to_test=2

# Read in inference functions
source(covid_get_path(inference_file))
source('set_up_covariates_and_data_11regions.R')
pars[['region_to_test']] = region_to_test
      pomp_object <- make_pomp_object_covid(
          n_regions = pars$n_regions,
          n_age_groups = n_age_groups,
          input_params = pars,
          delta_t = deltaT,
          contacts=pomp_contacts,
          population=population_list,
          beta_scales=beta_scales,
          beta_covar=beta_covariate_column,
          frac_underreported=fraction_underreported,
          dmeasure_Csnippet = dmeasure_snippet,
          rprocess_Csnippet = rprocess_snippet,
          rinit_Csnippet = rinit_snippet,
          data=data,
          fitstart=simstart,
          time_column='time',
          obsnames=observed_names,
          transformations=transformation,
          inherit_parameters=FALSE
          )

pf = pfilter(pomp_object, Np=5000)
print(pf@loglik)