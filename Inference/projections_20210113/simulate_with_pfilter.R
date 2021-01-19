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
jobid.true = as.numeric(args[2])
num_cores=as.numeric(args[3])
maxjobs=as.numeric(args[4])
region_to_test=arrayid

stopifnot(region_to_test %in% c(-1, 1:11))

pars$region_to_test = region_to_test
regtemp = region_to_test
pars$tmin = simstart

## read in plotting data

# Set fit pars to MLE
AIC = NA
if (use_mle == T){
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

  mle_frame = read.csv(sprintf('agg_mif_%s.csv', model_name)) %>%
    filter(region_to_test == regtemp) %>%
    filter(loglik == max(loglik)) %>%
    select(pars_to_fit, loglik) %>%
    unlist()
  for (p in pars_to_fit){
    pars[[p]] = mle_frame[[p]]
    print(paste(p, pars[[p]]))
  }

  AIC = 2 * length(pars_to_fit) - 2 * mle_frame[['loglik']]
}

pomp_object = simulate_pomp_covid(
  n_regions = pars$n_regions,
  n_age_groups = n_age_groups,
  nsim=1, 
  input_params = pars,
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
  obsnames = names(fitting_data %>% select(-time)),
  inference=T,
  dmeasure_Csnippet=dmeasure_snippet,
  data=fitting_data
)

print('Running mif')


library(doParallel)
cl <- makeCluster(num_cores)
t.0 = Sys.time()
registerDoParallel(cl)

foreach(i=1:maxjobs, .combine='rbind', 
  .packages=c('pomp', 'dplyr', 'reshape2', 'foreach')) %dopar%{
  pf = pfilter(pomp_object, Np=8000, filter.traj=T) # Runs a pfilter and saves a trajectory
  saveRDS(pf, sprintf('%s/pfilter_sim_reg_%s_%s_%s.rds', output_dir, region_to_test, jobid.true, i)) # Save the resulting pfilter for debugging/records

  # Sample a filtered trajectory of the desired states
  traj = filter.traj(pf)
  rm(pf)
  times = as.numeric(names(traj[1,1,])) # get the times

  d = data.frame(t(traj[,1,])) %>%
      select(grep(pattern=sprintf("_%s$", region_to_test), names(.))) %>%
      mutate(time = times,
             .id = i)
} -> trajectory

simtime = Sys.time() - t.0
print(simtime)
write.csv(trajectory %>% mutate(simtime = simtime), sprintf('%s/pfilter_trajectory_reg_%s.csv', output_dir, region_to_test), row.names=F)