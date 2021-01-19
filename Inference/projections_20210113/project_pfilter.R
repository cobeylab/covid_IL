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
# Read in inference functions
source('./input_file_specification_no_age.R')
source(covid_get_path(function_file))
source('./set_up_covariates_and_data_no_age.R')


# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1

# Reading input args
args = commandArgs(trailingOnly=TRUE)
arrayid = as.numeric(args[1])
num_cores=as.numeric(args[2])
starti=as.numeric(args[3])
endi=as.numeric(args[4])
region_to_test=arrayid
stopifnot(region_to_test %in% c(-1, 1:11))

pars$region_to_test = region_to_test
regtemp = region_to_test
pars$tmin = simstart
sim_times = data.frame(time=47:project_end)
pars$tmax = project_end

fitting_data = fitting_data %>%
  full_join(sim_times, by='time') %>%
  arrange(time)
finalobsnames = names(fitting_data %>% select(-time))
max_report_time = max(emr_report_covar$time)
covar_region =  emr_report_covar %>%
  select(time, grep(pattern=sprintf("_%s$", region_to_test), names(.))) %>%
  full_join(sim_times, by='time') %>%
  arrange(time)
covar_future = covar_region %>%
  filter(time > max_report_time) %>%
  mutate_at(vars(-time), ~0.99)
covar_region = covar_region %>%
  filter(time <= max_report_time) %>%
  bind_rows(covar_future) %>%
  arrange(time)


# Expand the slice dataframe
slice_df = read.csv(sprintf('slice_%s_%s.csv', model_name, region_to_test)) %>%
    filter(replicates > 0)
print(sprintf('slice_%s_%s.csv', model_name, region_to_test))
parstemp = pars
foreach (r=1:nrow(slice_df), .combine='rbind') %do%{
    n_sim_project = round(slice_df[r, 'replicates'])
    print(n_sim_project)
    slice_df[rep(r, n_sim_project), ]
} -> final_sample

# get pars to fit
reg_test = region_to_test
n_changepoints_reg = pars[[paste0('n_changepoints_', region_to_test)]]
n_hfr_changepoints_reg = pars[[paste0('n_HFR_changepoints_', region_to_test)]]
n_icu_changepoints_reg = pars[[paste0('n_ICU_changepoints_', region_to_test)]]
pars_to_fit = paste0(c('num_init_', 'phi_scale_', 
'inv_mu_0_', 'inv_mu_f_', 'inv_gamma_0_', 'inv_gamma_f_', 'inv_zeta_icu_0_', 'inv_zeta_icu_f_'), reg_test)
beta_par_names = paste0(paste0('beta_values_', 1:n_changepoints_reg), '_', reg_test)
hfr_par_names = paste0(paste0('HFR_values_', 1:n_hfr_changepoints_reg), '_', reg_test)
icu_par_names = paste0(paste0('ICU_values_', 1:n_icu_changepoints_reg), '_', reg_test)
pars_to_fit = c(pars_to_fit, beta_par_names, hfr_par_names, icu_par_names)

library(doParallel)
registerDoParallel(num_cores)

print(nrow(final_sample))

foreach(i=starti:endi, .combine='rbind'
  ) %dopar%{
  
  temppars = pars
  for (p in pars_to_fit){
    temppars[[p]] = as.numeric(final_sample[i, p])
  }
  pomp_object = simulate_pomp_covid(
        n_regions = temppars$n_regions,
        n_age_groups = n_age_groups,
        nsim=1, 
        input_params = temppars,
        delta_t = deltaT,
        population_list=population_list,
        death_reporting_covar=fraction_underreported,
        nonhosp_death_covar=nonhosp_deaths,
        beta_covar=beta_cov,
        beta_scale_covar=beta_scale_covar,
        emr_report_covar=emr_report_covar,
        rprocess_Csnippet = rprocess_snippet,
        rinit_Csnippet = rinit_snippet,
        rmeasure_Csnippet = rmeasure_snippet,
        obsnames = finalobsnames,
        inference=T,
        dmeasure_Csnippet=dmeasure_snippet,
        data=fitting_data
      )
  pf = pfilter(pomp_object, Np=5000, filter.traj=T) # Runs a pfilter and saves a trajectory
  saveRDS(pf, sprintf('%s/sim_pfilter_reg_%s_%s.csv', output_dir, region_to_test,i))
  # Sample a filtered trajectory of the desired states
  rm(pf)
  gc()
} -> trajectory

