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
library(doParallel)

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
args = commandArgs(trailingOnly=TRUE)
region_to_test=9#as.numeric(args[1])


stopifnot(region_to_test %in% c(-1, 1:11))
pars$region_to_test = region_to_test

n_changepoints_reg = unlist(pars[paste0('n_changepoints_', region_to_test)], use.names=F)
n_hfr_changepoints_reg = unlist(pars[paste0('n_HFR_changepoints_', region_to_test)], use.names=F)
n_icu_changepoints_reg = unlist(pars[paste0('n_ICU_changepoints_', region_to_test)], use.names=F)
partrans = parameter_trans(
    logit=c(
        paste0('num_init_', 1:pars$n_regions), 
        paste0('phi_scale_', 1:pars$n_regions),
        paste0('inv_mu_0_', 1:pars$n_regions),
        paste0('inv_mu_f_', 1:pars$n_regions),
        paste0('inv_gamma_0_', 1:pars$n_regions),
        paste0('inv_gamma_f_', 1:pars$n_regions),
        paste0('inv_zeta_icu_0_', 1:pars$n_regions),
        paste0('inv_zeta_icu_f_', 1:pars$n_regions),
        paste0('beta_values_', 1:n_changepoints_reg, '_', region_to_test),
        paste0('HFR_values_', 1:n_hfr_changepoints_reg, '_', region_to_test),
        paste0('ICU_values_', 1:n_icu_changepoints_reg, '_', region_to_test)
        )
)

time_pars = unlist(pars[paste0('changepoints_', 1:n_changepoints_reg, '_', region_to_test)], use.names=F)
par_names = paste0(paste0('beta_values_', 1:n_changepoints_reg), '_', region_to_test)
maxt=1000
foreach(i=1:(length(time_pars))) %do%{
    t1 = ifelse(i == 1, 0, time_pars[i - 1])
    t2 = ifelse(i == length(time_pars), maxt, time_pars[i + 1])
    parname = par_names[i]
    s = sprintf('%s=ifelse(time >= %s & time <= %s, %s, 0)', parname, t1, t2, sd_val)
    s
} -> strings
rw_beta = paste0(as.character(unlist(strings)), sep=', ', collapse='')

time_pars = unlist(pars[paste0('HFR_changepoint_values_', 1:n_hfr_changepoints_reg, '_', region_to_test)], use.names=F)
hfr_par_names = paste0(paste0('HFR_values_', 1:n_hfr_changepoints_reg), '_', region_to_test)
maxt=1000
foreach(i=1:(length(time_pars))) %do%{
    t1 = ifelse(i == 1, 0, time_pars[i - 1])
    t2 = ifelse(i == length(time_pars), maxt, time_pars[i + 1])
    parname = hfr_par_names[i]
    s = sprintf('%s=ifelse(time >= %s & time <= %s, %s, 0)', parname, t1, t2, sd_val)
    s
} -> strings
rw_hfr = paste0(as.character(unlist(strings)), sep=', ', collapse='')

time_pars = unlist(pars[paste0('ICU_changepoint_values_', 1:n_icu_changepoints_reg, '_', region_to_test)], use.names=F)
icu_par_names = paste0(paste0('ICU_values_', 1:n_icu_changepoints_reg), '_', region_to_test)
maxt=1000
foreach(i=1:(length(time_pars))) %do%{
    t1 = ifelse(i == 1, 0, time_pars[i - 1])
    t2 = ifelse(i == length(time_pars), maxt, time_pars[i + 1])
    parname = icu_par_names[i]
    s = sprintf('%s=ifelse(time >= %s & time <= %s, %s, 0)', parname, t1, t2, sd_val)
    s
} -> strings
rw_icu = paste0(as.character(unlist(strings)), sep=', ', collapse='')

pars_to_fit = paste0(c('num_init_', 'phi_scale_', 
  'inv_mu_0_', 'inv_mu_f_', 'inv_gamma_0_', 'inv_gamma_f_', 'inv_zeta_icu_0_', 'inv_zeta_icu_f_'), region_to_test)
rw_sd_str = do.call(sprintf, c(fmt="%s=ivp(0.05), %s=0.01, %s=ifelse(time<=147, 0.01, 0), %s=ifelse(time>=117, 0.01, 0), %s=ifelse(time<=147, 0.01, 0), %s=ifelse(time>=117, 0.01, 0), %s=ifelse(time<=147, 0.01, 0), %s=ifelse(time>=117, 0.01, 0)", as.list(pars_to_fit)))
rw_sd_str = paste(rw_beta, rw_hfr, rw_icu, rw_sd_str)
rw_sd_str = sprintf("rw.sd(%s)", rw_sd_str)
rw_sd = eval(parse(text=rw_sd_str))

print(fitting_data)

#region_fit_map = c(`10` = 11, `9` = 11, `1` = 2, `4`=5)

if(start_fit_at_mle == T){
  
  reg_test = region_to_test
  pars_to_fit = paste0(c('num_init_', 'phi_scale_', 
  'inv_mu_0_', 'inv_mu_f_', 'inv_gamma_0_', 'inv_gamma_f_', 'inv_zeta_icu_0_', 'inv_zeta_icu_f_'), reg_test)
  beta_par_names = paste0(paste0('beta_values_', 1:n_changepoints_reg), '_', reg_test)
  hfr_par_names = paste0(paste0('HFR_values_', 1:n_hfr_changepoints_reg), '_', reg_test)
  icu_par_names = paste0(paste0('ICU_values_', 1:n_icu_changepoints_reg), '_', reg_test)
  final_pars_to_fit = c(pars_to_fit, beta_par_names, hfr_par_names, icu_par_names)
  mle_frame = read.csv(mle_file)
  final_pars_to_fit = intersect(names(mle_frame), final_pars_to_fit)
  print(final_pars_to_fit)
  mle_frame = read.csv(mle_file) %>%
    filter(region_to_test == reg_test) %>%
    filter(loglik == max(loglik)) %>%
    select(final_pars_to_fit, loglik) %>%
    unlist()
  for(p in final_pars_to_fit){
    pars[[p]] = mle_frame[[p]]
  }
}
print(names(fitting_data %>% select(-time)))

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
      data=fitting_data,
      transformations=partrans
    )

pf = pfilter(pomp_object, Np=8000)
print(pf@loglik)