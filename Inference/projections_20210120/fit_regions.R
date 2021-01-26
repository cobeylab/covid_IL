## Run mif searches, can be implemented in parallel on the cluster 
## Sylvia Ranjeva, Phil Arevalo
## April 2020
## -----------------------------------------------------

## Load packages and specify files 
library(dplyr)
library(pomp)
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
# Read in inference functions
source('./input_file_specification_no_age.R')
source(covid_get_path(function_file))
source('./set_up_covariates_and_data_no_age.R')

## Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1
args = commandArgs(trailingOnly=TRUE)
arrayid = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
jobid = as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))
num_cores = as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
maxjobs=as.numeric(args[1])
region_to_test = arrayid

stopifnot(region_to_test %in% c(-1, 1:11))

fitting_data = fitting_data %>%
  arrange(time) %>%
  select(time, all_of(paste0(c('ObsHosp_', 'ObsDeaths_', 'ObsHospDeaths_', 'ObsICU_'), region_to_test)))
names(fitting_data) = sub(sprintf("_%s$", region_to_test), "_1", names(fitting_data))
finalobsnames = names(fitting_data %>% select(-time))

pars$region_to_test = region_to_test


# Set up transformations and rw.sd based on assumption that we're fitting regions separately
n_changepoints_reg = unlist(pars[paste0('n_changepoints_', region_to_test)], use.names=F)
n_hfr_changepoints_reg = unlist(pars[paste0('n_HFR_changepoints_', region_to_test)], use.names=F)
n_icu_changepoints_reg = unlist(pars[paste0('n_ICU_changepoints_', region_to_test)], use.names=F)
partrans = parameter_trans(
    logit=c(
        paste0('num_init_', 1), 
        paste0('phi_scale_', 1),
        paste0('inv_mu_0_', 1),
        paste0('inv_mu_f_', 1),
        paste0('inv_gamma_0_', 1),
        paste0('inv_gamma_f_', 1),
        paste0('inv_zeta_icu_0_', 1),
        paste0('inv_zeta_icu_f_', 1),
        paste0('beta_values_', 1:n_changepoints_reg, '_', 1),
        paste0('HFR_values_', 1:n_hfr_changepoints_reg, '_', 1),
        paste0('ICU_values_', 1:n_icu_changepoints_reg, '_', 1)
        )
)
rw_beta = get_rw_sd_str_time_var(time_pars = unlist(pars[paste0('changepoints_', 1:n_changepoints_reg, '_', region_to_test)], use.names=F), 
  par_names=paste0(paste0('beta_values_', 1:n_changepoints_reg), '_', 1), 
  sd_val)
rw_hfr =  get_rw_sd_str_time_var(time_pars = unlist(pars[paste0('HFR_changepoint_values_', 1:n_hfr_changepoints_reg, '_', region_to_test)], use.names=F), 
  par_names=paste0(paste0('HFR_values_', 1:n_hfr_changepoints_reg), '_', 1), 
  sd_val)
rw_icu = get_rw_sd_str_time_var(time_pars = unlist(pars[paste0('ICU_changepoint_values_', 1:n_icu_changepoints_reg, '_', region_to_test)], use.names=F), 
  par_names=paste0(paste0('ICU_values_', 1:n_icu_changepoints_reg), '_', 1), 
  sd_val)
pars_to_fit = paste0(c('num_init_', 'phi_scale_', 
  'inv_mu_0_', 'inv_mu_f_', 'inv_gamma_0_', 'inv_gamma_f_', 'inv_zeta_icu_0_', 'inv_zeta_icu_f_'), 1)
rw_sd_str = do.call(sprintf, c(fmt="%s=ivp(0.02), %s=0.02, %s=ifelse(time<=147, 0.02, 0), %s=ifelse(time>=117, 0.02, 0), %s=ifelse(time<=147, 0.02, 0), %s=ifelse(time>=117, 0.02, 0), %s=ifelse(time<=147, 0.02, 0), %s=ifelse(time>=117, 0.02, 0)", as.list(pars_to_fit)))
rw_sd_str = paste(rw_beta, rw_hfr, rw_icu, rw_sd_str)
rw_sd_str = sprintf("rw.sd(%s)", rw_sd_str)
rw_sd = eval(parse(text=rw_sd_str))


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
    filter(loglik == max(loglik, na.rm=T)) %>%
    select(all_of(final_pars_to_fit), loglik) %>%
    unlist()
  for(p in final_pars_to_fit){
    pars[[p]] = mle_frame[[p]]
  }
}


 pomp_object = simulate_pomp_covid(
      n_regions = pars$n_regions,
      n_age_groups = n_age_groups,
      nsim=1, 
      input_params = pars,
      delta_t = deltaT,
      population_list=population_list,
      death_reporting_covar=fraction_underreported,
      emr_report_covar=emr_report_covar,
      rprocess_Csnippet = rprocess_snippet,
      rinit_Csnippet = rinit_snippet,
      rmeasure_Csnippet = rmeasure_snippet,
      obsnames = finalobsnames,
      inference=T,
      dmeasure_Csnippet=dmeasure_snippet,
      data=fitting_data,
      transformations=partrans,
      global_Csnippet=global_snippet
    )

registerDoParallel(cores=num_cores)
foreach(i=1:maxjobs) %dopar%{
    mf = mif2(pomp_object,
    Nmif=n_mif,
    rw.sd = rw_sd,
    cooling.type="geometric",
    cooling.fraction.50=cooling_rate,
    Np=n_particles_mif)
  saveRDS(mf, sprintf("%s/reg_%s_mif_%s_%s.rds", output_dir, region_to_test, jobid, i))

  pf = pfilter(mf, Np=n_particles_pfilter)
  saveRDS(pf, sprintf("%s/reg_%s_pfilter_%s_%s.rds", output_dir, region_to_test, jobid, i))
} ->dummy