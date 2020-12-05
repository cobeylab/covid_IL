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
region_to_test=as.numeric(args[1])
stopifnot(region_to_test %in% c(-1, 1:11))

## read in plotting data
plotting_data = read.csv(covid_get_path(emr_data_file)) %>%
  filter(covid_region == region_to_test) %>%
  select(time, ObsHosp_1, ObsHospDeaths_1) %>%
  mutate(Date = as.Date('2020-01-14') + time) %>%
  arrange(Date) %>%
  mutate(Cumulative_Deaths=cumsum(if_else(is.na(ObsHospDeaths_1), 0, as.double(ObsHospDeaths_1)))) %>% 
  select(Date, Cumulative_Deaths, ObsHosp_1, ObsHospDeaths_1) %>%
  gather(Compartment, Cases, Cumulative_Deaths, ObsHosp_1, ObsHospDeaths_1) %>%
  mutate(Compartment = case_when((Compartment == 'ObsHospDeaths_1') ~ 'Reported hospital deaths',
                            (Compartment == 'ObsHosp_1') ~ 'Reported number of hospital beds occupied\nby COVID patients',
                            (Compartment == 'Cumulative_Deaths') ~ 'Reported cumulative hospital deaths'),
    Source = 'Observed') %>%
  select(Date, Compartment, Cases, Source) %>%
  mutate(Date = Date - 1)

ll_deaths = read.csv(covid_get_path(ll_file)) %>%
  filter(covid_region == region_to_test) %>%
  mutate(date = as.Date(date)) %>%
  rename(Date = date, Cases = ll_deaths) %>%
  mutate(Compartment = 'Reported total incident deaths') %>%
  mutate(Source='Linelist')

cumulative_ll = ll_deaths %>%
  arrange(Date) %>%
  mutate(Cases = cumsum(Cases),
    Compartment = 'Reported cumulative deaths')

cli_data = read.csv(covid_get_path(cli_file)) %>%
  mutate(Date=as.Date(date),
    Compartment='new_hospitalizations',
    Source='CLI') %>%
  rename(Cases = cli) %>%
  filter(covid_region == region_to_test) %>%
  select(Date, Cases, Compartment, Source)

plotting_data = plotting_data %>%
  bind_rows(ll_deaths) %>%
  bind_rows(cli_data) %>%
  bind_rows(cumulative_ll)

pars$region_to_test = region_to_test
regtemp = region_to_test
pars$tmin = simstart
pars$tmax = project_end

# Set fit pars to MLE
AIC = NA
if (use_mle == T){
  n_changepoints_reg = pars[[paste0('n_changepoints_', region_to_test)]]
  n_hfr_changepoints_reg = pars[[paste0('n_HFR_changepoints_', region_to_test)]]
  pars_to_fit = paste0(c('num_init_', 'phi_scale_', 'inv_mu_0_', 'inv_mu_f_', 'inv_gamma_0_', 'inv_gamma_f_'), region_to_test)
  beta_par_names = paste0(paste0('beta_values_', 1:n_changepoints_reg), '_', region_to_test)
  hfr_par_names = paste0(paste0('HFR_values_', 1:n_hfr_changepoints_reg), '_', region_to_test)
  pars_to_fit = c(pars_to_fit, beta_par_names, hfr_par_names)
  print(pars_to_fit)
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


## Simulate
sim_full = simulate_pomp_covid(
  n_regions = pars$n_regions,
  n_age_groups = n_age_groups,
  nsim=n_sim, 
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
  obsnames = c(paste0('ObsHosp_', 1:11), paste0('ObsHospDeaths_',1:11), paste0('ObsDeaths_', 1:11))
)



## Process output

if (region_to_test == -1){
  region_criteria = 1:11
} else{
  region_criteria = c(region_to_test)
}

pops = population_list %>%
  group_by(covid_region) %>%
  summarize(pop=sum(POPULATION)) %>%
  ungroup() 
popfinal = pops %>% select(pop) %>% unlist()
names(popfinal) = pops$covid_region

alloutputs = process_pomp_covid_output(sim_full, agg_regions=F)

print("Preparing for plotting")
plotout_noage = alloutputs$plotting_output %>% ungroup() %>%    
    mutate(Date=as.Date('2020-01-14') + Time) %>%
    filter(Compartment %in% c('ObsHospDeaths', 'ObsDeaths', "ObsHosp", 'Inc', "Total infectious",  "R", 'new_hospitalizations', 
      'transmissionRate', 'HFRtrack','IHRtrack','phitrack'),
      Region %in% region_criteria) %>%
    mutate(Compartment = case_when((Compartment == 'ObsHospDeaths') ~ 'Reported hospital deaths',
                            (Compartment == 'ObsHosp') ~ 'Reported number of hospital beds occupied\nby COVID patients',
                            (Compartment == 'new_hospitalizations') ~ 'new_hospitalizations',
                            (Compartment == 'Inc') ~ "Daily new infections",
                            (Compartment == 'Total infectious') ~ "Total infectious",
                            (Compartment == 'ObsDeaths') ~ 'Reported total incident deaths',
                            (Compartment == "R") ~ 'Fraction recovered',
                            (Compartment == 'transmissionRate') ~ 'Transmission rate',
                            (Compartment == 'phitrack') ~ 'Fraction of deaths outside the hospital',
                            (Compartment == 'HFRtrack') ~ 'HFR',
                            (Compartment == 'IHRtrack') ~ 'IHR',
                            T ~ Compartment),
      Source = "Simulation") %>%
    group_by(Region, Date, Compartment, SimID) %>%
    dplyr::summarize(Cases=sum(Cases)) %>%
    mutate(Source = 'Simulation') %>%
    ungroup()

plot_order = c('Reported hospital deaths','Reported cumulative hospital deaths', 'Reported number of hospital beds occupied\nby COVID patients',
  'Reported total incident deaths', 'Reported cumulative deaths', 'new_hospitalizations', 
  'Daily new infections', 'Total infectious', 'Fraction recovered', 
  'HFR','IHR','Fraction of deaths outside the hospital',
  'IFR', 'Transmission rate', 'R(t)',
  'Hospital duration (recover)', 'Hospital duration (death)')

compartment_order = c('Observed', 'Simulation', 'Linelist', 'CLI')

#head(alloutputs$plotting_output %>% filter(Compartment == 'IH', Region=='1'))
#h = sum(plotout_noage %>% filter(Compartment == 'new_hospitalizations') %>% select(Cases) %>% unlist())
#infe = sum(plotout_noage %>% filter(Compartment == 'Incidence') %>% select(Cases) %>% unlist(), na.rm=T)
#hdeath = sum(plotout_noage %>% filter(Compartment == 'Hospital deaths') %>% select(Cases) %>% unlist(), na.rm=T)
#print(h/infe)
#print(hdeath/h)

cumulative_hosp = plotout_noage %>%  
    filter(Compartment=="Reported hospital deaths") %>%
    group_by(SimID) %>%
    arrange(Date) %>%
    mutate(Cases=cumsum(Cases)) %>%
    ungroup() %>%
    mutate(Compartment='Reported cumulative hospital deaths')

cumulative_all_deaths = plotout_noage %>%  
    filter(Compartment=="Reported total incident deaths") %>%
    group_by(SimID) %>%
    arrange(Date) %>%
    mutate(Cases=cumsum(Cases)) %>%
    ungroup() %>%
    mutate(Compartment='Reported cumulative deaths')

D= alloutputs$plotting_output %>% ungroup() %>%  
    filter(Compartment=="D", Region %in% region_criteria) %>%
    group_by(SimID, Time) %>%
    summarize(Dead=sum(Cases)) %>%
    ungroup() %>%
    mutate(Time = Time - 19)

R= alloutputs$plotting_output %>% ungroup() %>%  
    filter(Compartment=="R", Region %in% region_criteria) %>%
    group_by(SimID, Time) %>%
    summarize(Recovered=sum(Cases)) %>%
    ungroup() %>%
    mutate(Time = Time - 9)

IFRdf = full_join(D, R, by=c('SimID', 'Time')) %>%
    drop_na() %>%
    mutate(IFR = Dead /(Dead+Recovered)) %>%
    select(SimID, Time, IFR) %>%
    gather(Compartment, Cases, IFR) %>%
    mutate(Date=as.Date('2020-01-14') + Time, 
    Source='Simulation') %>%
    drop_na() %>%
    filter(Date>= as.Date('2020-03-20'))

plotout_noage = bind_rows(plotout_noage, IFRdf) %>%
  bind_rows(cumulative_hosp) %>%
  bind_rows(cumulative_all_deaths)

frac_recovered = plotout_noage %>%
  filter(Compartment == 'Fraction recovered') %>%
  mutate(Cases = Cases / popfinal[Region])

beta_init = unlogit(pars[[paste0('beta_values_1_', region_to_test)]], 0.9, 0.4)
HFR_init = unlogit(pars[[paste0('HFR_values_1_', region_to_test)]], pars$HFR_max, pars$HFR_min)
phi = pars[[paste0('phi_scale_', region_to_test)]] * HFR_init
IHR_init = (pars$IFR_constraint - phi) / (HFR_init - phi)
paras = list(frac_nonhosp_deaths = (phi * (1-IHR_init)) / (phi * (1-IHR_init) + IHR_init * HFR_init),
                     sigma = 1/3.5,
                     zeta_s = 1/1.5,
                     mu_m = 1/15,
                     gamma_m = 1/3.5,
                     zeta_h = 1/4.2)
r0pop = read.csv(covid_get_path('Data/covid_region_populations.csv')) %>% 
    summarize(POPULATION=sum(POPULATION))

source('R0_functions.R')
R0_val = get_R0(region_cons=region_to_test, 
       beta_value=beta_init,
       pop=r0pop,
       paras=paras,
       t=47,
       IHR =IHR_init,
       HFR=HFR_init)

Rt_scaling_constant = 5.05679
plotout_noage = bind_rows(plotout_noage %>% filter(Compartment != "Fraction recovered"), frac_recovered)
Rt = plotout_noage %>%
  filter(Compartment == 'Transmission rate' | Compartment == 'Fraction recovered') %>%
  spread(Compartment, Cases) %>%
  mutate(Rt = Rt_scaling_constant * `Transmission rate` * (1 - `Fraction recovered`)) %>%
  select(-`Transmission rate`, -`Fraction recovered`) %>%
  rename(Cases=Rt) %>%
  mutate(Compartment='R(t)')
plotout_noage = bind_rows(plotout_noage, Rt)

## add in gamma and mu
g0 = unlogit(pars[[paste0('inv_gamma_0_', region_to_test)]], pars$gamma_max, pars$gamma_min)
gf = unlogit(pars[[paste0('inv_gamma_f_', region_to_test)]], pars$gamma_max, pars$gamma_min)

m0 = unlogit(pars[[paste0('inv_mu_0_', region_to_test)]], pars$mu_max, pars$mu_min)
mf = unlogit(pars[[paste0('inv_mu_f_', region_to_test)]], pars$mu_max, pars$mu_min)

times = sort(unique(as.numeric(plotout_noage$Date - as.Date('2020-01-14'))))

g = mapply(FUN=get_par_value, times, p0=g0, pf=gf, t0=pars$hosp_t_min, tf=pars$hosp_t_max)
m = mapply(FUN=get_par_value, times, p0=m0, pf=mf, t0=pars$hosp_t_min, tf=pars$hosp_t_max)

gamma_df = data.frame(Date = as.Date('2020-01-14') + times, Compartment='Hospital duration (recover)', Cases=g, Source='Simulation')
mu_df = data.frame(Date = as.Date('2020-01-14') + times, Compartment='Hospital duration (death)', Cases=m, Source='Simulation')


#pardf = data.frame(Date = unique(Date(plotout_noage))) %>%
 # mutate(Time = Date - as.Date('2020-01-14'),
#    mu_h = )

plotout_noage = bind_rows(plotout_noage, gamma_df) %>% bind_rows(mu_df) %>%
  filter(Compartment %in% c('Reported number of hospital beds occupied\nby COVID patients', 'Fraction recovered')) %>%
  filter(Date >= as.Date('2020-06-01'))

plotting_data = plotting_data %>% 
  filter(Compartment %in% c('Reported number of hospital beds occupied\nby COVID patients', 'Fraction recovered')) %>%
  filter(Date >= as.Date('2020-06-01'))

print("Plotting")
ggplot(plotout_noage, aes(x=Date, y=Cases, color=Source, fill=Source)) +
   stat_summary(fun=median, geom="line", size=0.75) +
   stat_summary(fun.min=function(z){quantile(z,0.025)}, fun.max=function(z){quantile(z,0.975)}, geom="ribbon", alpha=0.3, color=NA) +
   geom_point(plotting_data, mapping=aes(x=Date, y=Cases, fill=Source), size=0.9, alpha=0.7, pch=21) +
   facet_wrap(~factor(Compartment, levels=plot_order), scales='free_y', ncol=3) +
   theme_bw() +
   ylab('') +
   scale_x_date(date_breaks = "months", date_labels="%b") +
   ggtitle(sprintf("Chicago", region_to_test, R0_val, round(popfinal[[region_to_test]] * pars[[paste0('num_init_', region_to_test)]]), AIC)) +
   scale_color_brewer(palette='Set1') +
   scale_fill_brewer(palette='Set1') +
   labs(fill='', color='') +
   geom_hline(data.frame(Compartment = 'Reported number of hospital beds occupied\nby COVID patients', yint=2442), mapping=aes(yintercept=yint))
ggsave(sprintf('%s_%s.png', model_name, region_to_test), width=8, height=3)


#print("Writing output")
#write.csv(alloutputs$plotting_output, 'ploutout.csv', row.names=F)