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
source(covid_get_path(icu_function_file))
source('./set_up_covariates_and_data_no_age.R')


# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1

# Reading input args
args = commandArgs(trailingOnly=TRUE)
region_to_test=as.numeric(args[1])
stopifnot(region_to_test %in% c(-1, 1:11))

## read in plotting data
plotting_data = read.csv(covid_get_path(emr_data_file)) %>%
  filter(covid_region == region_to_test) %>%
  select(time, ObsICU_1, ObsHospDeaths_1) %>%
  mutate(Date = as.Date('2020-01-14') + time) %>%
  arrange(Date) %>%
  select(Date, ObsICU_1) %>%
  gather(Compartment, Cases, ObsICU_1) %>%
  mutate(Compartment = 'ICU beds occupied',
    Source = 'EMResource') %>%
  select(Date, Compartment, Cases, Source) %>%
  mutate(Date = Date)
head(plotting_data)
pars$region_to_test = region_to_test
regtemp = region_to_test
pars$tmin = simstart

# Set fit pars to MLE
AIC = NA
if (use_mle == T){
  n_changepoints_reg = pars[[paste0('n_changepoints_', region_to_test)]]
  n_hfr_changepoints_reg = pars[[paste0('n_HFR_changepoints_', region_to_test)]]
  pars_to_fit = paste0(c('num_init_', 'phi_scale_', 'inv_mu_0_', 'inv_mu_f_', 'inv_gamma_0_', 'inv_gamma_f_'), region_to_test)
  beta_par_names = paste0(paste0('beta_values_', 1:n_changepoints_reg), '_', region_to_test)
  hfr_par_names = paste0(paste0('HFR_values_', 1:n_hfr_changepoints_reg), '_', region_to_test)
  pars_to_fit = c(pars_to_fit, beta_par_names, hfr_par_names)
  mle_frame = read.csv(sprintf('agg_mif_%s.csv', model_name)) %>%
    filter(region_to_test == regtemp) %>%
    filter(loglik == max(loglik)) %>%
    select(pars_to_fit, loglik) %>%
    unlist()


  for (p in pars_to_fit){
    pars[[p]] = mle_frame[[p]]
  }

  AIC = 2 * length(pars_to_fit) - 2 * mle_frame[['loglik']]
}


inchosp_full = read.csv('projection_outputs_beta_hfr_mu_gamma_constphi_cdc_outputs/projections_beta_hfr_mu_gamma_constphi_cdc_11.csv') %>%
  filter(Compartment == 'new_hospitalizations') %>%
  mutate(time=as.numeric(as.Date(date) - as.Date('2020-01-14')),
    new_hospitalizations=Cases) %>%
  select(time, new_hospitalizations, SimID)

## Simulate

foreach(s=1:100, .combine='rbind') %do%{

  inchosp = inchosp_full %>% filter(SimID==s)
  sim_full = simulate_pomp_covid(
    n_regions = pars$n_regions,
    n_age_groups = n_age_groups,
    nsim=1, 
    input_params = pars,
    delta_t = deltaT,
    emr_report_covar=emr_report_covar,
    icu_covars=icu_covar,
    inc_hosp=inchosp,
    rprocess_Csnippet = rprocess_icu_snippet,
    rinit_Csnippet = rinit_icu_snippet
  ) 
    outputs = process_pomp_covid_output(sim_full$raw_simulation_output, agg_regions=F)
    outputs$plotting_output %>% ungroup() %>% mutate(SimID=s)
} -> alloutputs


head(alloutputs)

plotout_noage = alloutputs %>% ungroup() %>%    
    mutate(Date=as.Date('2020-01-14') + Time) %>%
    filter(Compartment %in% c('IC'),
      Region == region_to_test) %>%
    mutate(
      Source = 'Simulation',
      Compartment = 'ICU beds occupied') %>%
    group_by(Region, Date, Compartment, SimID, Source) %>%
    dplyr::summarize(Cases=sum(Cases)) %>%
    ungroup()



print("Plotting")
ggplot(plotout_noage, aes(x=Date, y=Cases, color=Source, fill=Source)) +
   stat_summary(fun=function(z){quantile(z,0.5,type=3)}, geom="line", size=0.75) +
   stat_summary(fun.min=function(z){quantile(z,0.025, type=3)}, fun.max=function(z){quantile(z,0.975, type=3)}, geom="ribbon", alpha=0.3, color=NA) +
   geom_point(plotting_data, mapping=aes(x=Date, y=Cases, fill=Source), size=0.9, alpha=0.7, pch=21) +
   theme_bw() +
   ylab('') +
   scale_x_date(date_breaks = "months", date_labels="%b") +
   scale_color_brewer(palette='Set1') +
   scale_fill_brewer(palette='Set1')
ggsave(sprintf('%s_icu_%s.png', model_name, region_to_test), width=14, height=10)


#print("Writing output")
#write.csv(alloutputs$plotting_output, 'ploutout.csv', row.names=F)