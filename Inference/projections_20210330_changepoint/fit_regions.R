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
region_to_test = arrayid
maxjobs=as.numeric(args[1])
ifr_constraint=as.numeric(args[2])

stopifnot(region_to_test %in% c(-1, 1:11))

fitting_data = fitting_data %>%
  arrange(time) %>%
  select(time, all_of(paste0(c('ObsHosp_', 'ObsDeaths_', 'ObsICU_'), region_to_test)))
names(fitting_data) = sub(sprintf("_%s$", region_to_test), "_1", names(fitting_data))
finalobsnames = names(fitting_data %>% select(-time))

pars$region_to_test = region_to_test
print(sprintf("%s/reg_%s_pfilter_*.rds", previous_inference, region_to_test))
if(start_fit_at_mle == T){
  files = Sys.glob(sprintf("%s/reg_%s_pfilter_*.rds", previous_inference, region_to_test))
  print(sprintf("%s/reg_%s_pfilter_*.rds", previous_inference, region_to_test))
  best = files[1]
  loglik.max = -1e10
  for (f in files){
    pf = readRDS(f)
    if (pf@loglik > loglik.max & pf@params[['IFR_constraint']]==ifr_constraint){
      best = f
      loglik.max = pf@loglik
    }
  }
  pf = readRDS(best)
  parnames=names(pf@params)
  newpars = as.numeric(pf@params)
  names(newpars) = parnames

#  n_transmission_changepoints = newpars[['n_changepoints_1']] 

  #beta_indices = grep('^beta_values*', names(newpars))
  #changepoint_indices = grep('^changepoints*', names(newpars))

  #beta_values = newpars[beta_indices]
  #beta_values[[sprintf('beta_values_%s_1', n_transmission_changepoints + 1)]] =  beta_values[[sprintf('beta_values_%s_1', n_transmission_changepoints)]]
  #beta_values[[sprintf('beta_values_%s_1', n_transmission_changepoints + 2)]] =  beta_values[[sprintf('beta_values_%s_1', n_transmission_changepoints)]]

  #changepoints = newpars[changepoint_indices]
  #changepoints[[sprintf('changepoints_%s_1', n_transmission_changepoints + 1)]] =  changepoints[[sprintf('changepoints_%s_1', n_transmission_changepoints)]] + 15
  #changepoints[[sprintf('changepoints_%s_1', n_transmission_changepoints + 2)]] =  changepoints[[sprintf('changepoints_%s_1', n_transmission_changepoints)]] + 30
  
  #newpars[['n_changepoints_1']] = n_transmission_changepoints + 2

  #newpars = c(newpars[-c(beta_indices, changepoint_indices)], beta_values, changepoints)
}


newpars[['inv_mu_h_1']] = 9.598

# Set up transformations and rw.sd based on assumption that we're fitting regions separately
n_changepoints_reg = unlist(newpars[paste0('n_changepoints_', 1)], use.names=F)
n_hfr_changepoints_reg = unlist(pars[paste0('n_HFR_changepoints_', region_to_test)], use.names=F)
n_icu_changepoints_reg = unlist(pars[paste0('n_ICU_changepoints_', region_to_test)], use.names=F)
partrans = parameter_trans(
    logit=c(
        paste0('num_init_', 1), 
        paste0('beta_values_', 1:n_changepoints_reg, '_', 1),
        paste0('HFR_values_', 1:n_hfr_changepoints_reg, '_', 1),
        paste0('ICU_values_', 1:n_icu_changepoints_reg, '_', 1)#,
        #'init_freq'
        )
)
rw_beta = get_rw_sd_str_time_var(time_pars = unlist(newpars[paste0('changepoints_', 1:n_changepoints_reg, '_', 1)], use.names=F), 
  par_names=paste0(paste0('beta_values_', 1:n_changepoints_reg), '_', 1), 
  sd_val)
rw_hfr =  get_rw_sd_str_time_var(time_pars = unlist(pars[paste0('HFR_changepoint_values_', 1:n_hfr_changepoints_reg, '_', region_to_test)], use.names=F), 
  par_names=paste0(paste0('HFR_values_', 1:n_hfr_changepoints_reg), '_', 1), 
  sd_val)
rw_icu = get_rw_sd_str_time_var(time_pars = unlist(pars[paste0('ICU_changepoint_values_', 1:n_icu_changepoints_reg, '_', region_to_test)], use.names=F), 
  par_names=paste0(paste0('ICU_values_', 1:n_icu_changepoints_reg), '_', 1), 
  sd_val)

init_params = paste0(c('num_init_'), 1,'=ivp(', sd_val_init, ')', collapse=',')
extended_params = paste0(c('inv_gamma_h_intercept_', 'inv_zeta_icu_intercept_'), 1, '=', sd_val, collapse=',')
rw_sd_str = paste(init_params, extended_params, sep=',')


rw_sd_str = paste(rw_beta, rw_hfr, rw_icu, rw_sd_str)
rw_sd_str = sprintf("rw.sd(%s)", rw_sd_str)
rw_sd = eval(parse(text=rw_sd_str))


registerDoParallel(cores=num_cores)

allpars = expand.grid(sc=c(1.3, 1.5), 
  scd=c(1.3, 1.7), 
  all_init_freqs=seq(0.01, 0.99, length.out=4))

foreach(i=0:((nrow(allpars) * 3) - 1)) %dopar%{
  newpars[['scalefactor']] = allpars[i %% nrow(allpars) + 1, 'sc']
  newpars[['scalefactor_death']] = allpars[i %% nrow(allpars) + 1, 'scd']
  newpars[['init_freq']] = allpars[i %% nrow(allpars) + 1, 'all_init_freqs']
  pomp_object = simulate_pomp_covid(
      n_regions = pars$n_regions,
      n_age_groups = n_age_groups,
      nsim=1, 
      input_params = pars,
      delta_t = deltaT,
      population_list=population_list,
      region_covars=region_covars,
      shared_covars=shared_covars,
      rprocess_Csnippet = rprocess_scenario_snippet,
      rinit_Csnippet = rinit_snippet,
      rmeasure_Csnippet = rmeasure_snippet,
      obsnames = finalobsnames,
      inference=T,
      dmeasure_Csnippet=dmeasure_snippet,
      data=fitting_data,
      transformations=partrans,
      global_Csnippet=global_snippet,
      use.mle.params = start_fit_at_mle,
      mle.params=newpars
    )



  mf = mif2(pomp_object,
    Nmif=n_mif,
    rw.sd = rw_sd,
    cooling.type="geometric",
    cooling.fraction.50=cooling_rate,
    Np=n_particles_mif)
  saveRDS(mf, sprintf("%s/reg_%s_mif_%s_%s.rds", output_dir, region_to_test, jobid, i))

  for (k in 1:n_reps_pfilter){
      pf = pfilter(mf, Np=n_particles_pfilter)
      saveRDS(pf, sprintf("%s/reg_%s_pfilter_%s_%s_%s.rds", output_dir, region_to_test, jobid, i, k))

      #if ( i <= (nrow(allpars)  - 1)) {
      #    pf = pfilter(pomp_object, Np=n_particles_pfilter)
      #    saveRDS(pf, sprintf("%s/reg_%s_pfilter_%s_%s_%s_%s.rds", output_dir, region_to_test, jobid, i, 'orig', k))
      #}
  }

} ->dummy