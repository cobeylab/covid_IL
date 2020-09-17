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

root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

source('./input_file_specification.R')
source(covid_get_path(inference_file))
output_dir = './'
deltaT = 0.1
previous_mle_file = './orig.final_mle_pars.csv'
current_mle_file = './final_mle_pars.csv'


default_par_file = previous_mle_file
source('set_up_covariates_and_data.R')
pomp_object_previous <- make_pomp_object_covid(
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

default_par_file = current_mle_file
print('Initializing parameters')
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)
pars$n_regions = length(unique(population_list$region))

## Add in initial age distribution
age_dist_frame = read.csv(covid_get_path(age_dist_file)) 
temp_age_dist = age_dist_frame %>% 
  filter(b_elderly == pars$b_elderly) %>% select(-b_elderly)
age_dist = as.numeric(temp_age_dist$value)
names(age_dist) = temp_age_dist$param_name
pars = c(age_dist, pars)
pars$tmax = simend

pomp_object_current <- make_pomp_object_covid(
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

cl <- makeCluster(6)
registerDoParallel(cl)

Np = 5000
foreach(region_to_test = 1:6, .combine='rbind',
    .packages=c('pomp', 'dplyr', 'reshape2')) %dopar%{
    select <- dplyr::select
    rename <- dplyr::rename
    summarize <- dplyr::summarise
    contains <- dplyr::contains

    pomp_object_current@params[['region_to_test']] = region_to_test
    pomp_object_previous@params[['region_to_test']] = region_to_test
    pf_current = pfilter(pomp_object_current, Np = Np)
    pf_previous = pfilter(pomp_object_previous, Np = Np)
    c(region_to_test, pf_current@loglik, pf_previous@loglik)
} -> result
stopCluster(cl)
result = data.frame(result)
print(result)
names(result) = c('Region', 'Current MLE', 'Previous MLE')
write.csv(result, "mle_comparison.csv", row.names=F)
print(data.frame(result))