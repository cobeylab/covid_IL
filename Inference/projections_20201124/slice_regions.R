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
arrayid = as.numeric(args[1])
jobid.true = as.numeric(args[2])
num_cores=as.numeric(args[3])
maxjobs=as.numeric(args[4])
region_to_test=as.numeric(args[5])


stopifnot(region_to_test %in% c(-1, 1:11))
pars$region_to_test = region_to_test


loopstart = 1
loopend = maxjobs
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Set parameters to MLE
n_changepoints_reg = pars[[paste0('n_changepoints_', region_to_test)]]
n_hfr_changepoints_reg = pars[[paste0('n_HFR_changepoints_', region_to_test)]]
pars_to_fit = paste0(c('num_init_', 'phi_scale_', 'inv_mu_0_', 'inv_mu_f_', 'inv_gamma_0_', 'inv_gamma_f_'), region_to_test)
beta_par_names = paste0(paste0('beta_values_', 1:n_changepoints_reg), '_', region_to_test)
hfr_par_names = paste0(paste0('HFR_values_', 1:n_hfr_changepoints_reg), '_', region_to_test)
pars_to_fit = c(pars_to_fit, beta_par_names, hfr_par_names)

reg_test = region_to_test
mle_frame = read.csv(sprintf('agg_mif_%s.csv', model_name)) %>%
  filter(region_to_test == reg_test) %>%
  filter(loglik == max(loglik)) %>%
  select(pars_to_fit, loglik) %>%
  unlist()
for (p in pars_to_fit){
  pars[[p]] = mle_frame[[p]]
}

slice.par = sprintf('beta_values_%s_%s', n_changepoints_reg, region_to_test)
mle.par = pars[[slice.par]]

slice.values = seq(mle.par * 0.9, mle.par * 1.1, length.out=28)


foreach(jobid=loopstart:loopend,
  .packages=c('pomp', 'dplyr', 'reshape2'),
  .combine='rbind') %dopar% {
      stopifnot(jobid <= length(slice.values))

      pars[[slice.par]] = slice.values[jobid]

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

    # Time to run one mif with 3000 particles ~0.5 min
    # Time to run one mif with 5000 particles ~1.2 min

    # Time to run 5 replicate pfilters with 5000 particles ~ 5 mins
    pf = run_pfilter_and_output_result(
      n_reps_pfilter=n_reps_pfilter, 
      n_particles_pfilter=n_particles_pfilter, 
      jobid=paste0('slice_', jobid.true, '_', arrayid),
      arrayid=jobid,
      output_dir=output_dir,
      pomp_object=pomp_object,
      t_mif=NA,
      return_df = T)
   
  } -> slice.liks


#####
library(matrixStats)
mle_frame = read.csv(sprintf('agg_mif_%s.csv', model_name)) %>%
  filter(region_to_test == reg_test) %>%
  filter(loglik == max(loglik)) %>%
  bind_rows(slice.liks) %>%
  mutate(probs = exp(loglik - logSumExp(loglik)),
    replicates = rmultinom(1, n_sim, probs)) %>%
  arrange(loglik)

write.csv(mle_frame, sprintf("./slice_%s_%s.csv", model_name, region_to_test))
