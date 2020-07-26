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
library(doParallel)

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

# Set mif search parameters
n_mif = 50
n_particles_mif = 3000
cooling_rate = 0.95
num_points = 100

n_reps_pfilter = 3
n_particles_pfilter = 5000

root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

source('./input_file_specification.R')
default_par_file = './default_parameter_values.csv'

# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1

# Reading input args
args = commandArgs(trailingOnly=TRUE)
jobid_master = 1
arrayid = 1
num_cores=1
maxjobs=1
region_to_test=1
output_dir='./test'


# Read in inference functions
source(covid_get_path(inference_file))
source('set_up_covariates_and_data.R')


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