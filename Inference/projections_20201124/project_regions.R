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
library(purrr)

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
                            (Compartment == 'ObsHosp_1') ~ 'Reported census hospital beds occupied',
                            (Compartment == 'Cumulative_Deaths') ~ 'Reported cumulative hospital deaths'),
    Source = 'EMResource') %>%
  select(Date, Compartment, Cases, Source) %>%
  mutate(Date = Date - 1)

ll_deaths = read.csv(covid_get_path(ll_file)) %>%
  filter(covid_region == region_to_test) %>%
  mutate(date = as.Date(date)) %>%
  rename(Date = date, Cases = ll_deaths) %>%
  mutate(Compartment = 'Total incident deaths') %>%
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

sim_0 = 0

slice_df = read.csv(sprintf('slice_%s_%s.csv', model_name, region_to_test)) %>%
  filter(replicates > 0)

foreach (r=1:nrow(slice_df), .combine='rbind') %do%{
  n_sim_project = slice_df[r, 'replicates']
  for (p in pars_to_fit){
    pars[[p]] = slice_df[r, p]
  }

  ## Simulate
  sim_temp = simulate_pomp_covid(
    n_regions = pars$n_regions,
    n_age_groups = n_age_groups,
    nsim=n_sim_project, 
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
  sim_temp$raw_simulation_output = sim_temp$raw_simulation_output %>%
    mutate(.id = as.numeric(.id), .id = .id + sim_0)
  sim_0 = sim_0 + n_sim_project
  sim_temp$raw_simulation_output
} -> sim_full

print(length(unique(sim_full$.id)))

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

print("Writing output")
## Output civis format
projections = alloutputs$plotting_output %>%  
    filter(Region == region_to_test) %>%
    mutate(date=as.Date('2020-01-14') + Time) %>%
    select(-Time)
write.csv(projections, sprintf("%s/projections_%s_%s.csv", human_output_dir, model_name, region_to_test), row.names=F)

compartments = c(`Total infectious`='cases_',
  Inc='cases_new_',
  new_deaths='deaths_',
  ObsDeaths='deaths_det_',
  IH='hosp_bed_',
  R='recovered_')
foreach (i=1:length(compartments)) %do%{
    comp = names(compartments[i])
    comp_name = compartments[i]
    
    cases = projections %>%
        filter(date >= as.Date('2020-03-13')) %>%
        filter(Compartment == comp) %>%
        group_by(date) %>%
        dplyr::summarize(med=quantile(Cases, 0.5, type=3),
                         lower=quantile(Cases, 0.025, type=3),
                         upper=quantile(Cases, 0.975, type=3)) %>%
        mutate(geography_modeled = sprintf('covidregion_%s', region_to_test),
               scenario_name = 'baseline') %>%
        ungroup() %>%
        select(date, geography_modeled, scenario_name, med, lower, upper)
    names(cases) = c('date','geography_modeled','scenario_name',paste0(comp_name, c('median','lower','upper')))
    cases
} -> tempframe
finaldf = tempframe %>% reduce(full_join, by = c("date", "scenario_name", "geography_modeled"))
write.csv(finaldf, sprintf('%s/uchicago_covidregion_%s_%s.csv', human_output_dir, region_to_test, Sys.Date()-1), row.names=F)

## Output hospital capacity report
hb_template = read.csv(covid_get_path(hosp_capacity_file)) %>%
  filter(geography_modeled == paste0('covidregion_', region_to_test)) %>%
  filter(resource_type=='hb_availforcovid')
nsim = length(unique(projections$SimID))
for (r in 1:nrow(hb_template)){
    beds = as.numeric(hb_template[r, 'avg_resource_available'])
    d = as.Date(hb_template[r, 'date_window_upper_bound'])
    temp = projections %>% 
        filter(date <= d, d>=Sys.Date()-1) %>%
        filter(Compartment == 'IH') %>%
        group_by(SimID) %>%
        dplyr::summarize(exceeded = any(Cases >= beds)) %>%
        select(exceeded) %>%
        unlist(use.names=F)
    hb_template[r, 'percent_of_simulations_that_exceed'] = sum(temp)/nsim * 100
}
write.csv(hb_template %>% mutate(scenario_name='baseline'), sprintf('%s/uchicago_hospitaloverflow_covidregion_%s_%s.csv', human_output_dir, region_to_test, Sys.Date()-1), row.names=F)
rm(projections)


print("Preparing for plotting")
plotout_noage = alloutputs$plotting_output %>% ungroup() %>%    
    mutate(Date=as.Date('2020-01-14') + Time) %>%
    filter(Compartment %in% c('ObsHospDeaths', 'ObsDeaths', "ObsHosp", 'Inc', "Total infectious",  "R", 'new_hospitalizations', 
      'transmissionRate', 'HFRtrack','IHRtrack','phitrack', 'new_deaths'),
      Region %in% region_criteria) %>%
    mutate(
      Source = ifelse(Compartment %in% c('ObsHospDeaths', 'ObsDeaths', 'ObsHosp'), 'Simulation (reported)', 'Simulation'),
      Compartment = case_when((Compartment == 'ObsHospDeaths') ~ 'Reported hospital deaths',
                            (Compartment == 'ObsHosp') ~ 'Reported census hospital beds occupied',
                            (Compartment == 'new_hospitalizations') ~ 'new_hospitalizations',
                            (Compartment == 'Inc') ~ "Daily new infections",
                            (Compartment == 'Total infectious') ~ "Total infectious",
                            (Compartment == 'ObsDeaths') ~ 'Total incident deaths',
                            (Compartment == 'new_deaths') ~ 'Total incident deaths',
                            (Compartment == "R") ~ 'Fraction recovered',
                            (Compartment == 'transmissionRate') ~ 'Transmission rate',
                            (Compartment == 'phitrack') ~ 'Fraction of deaths outside the hospital',
                            (Compartment == 'HFRtrack') ~ 'HFR',
                            (Compartment == 'IHRtrack') ~ 'IHR',
                            T ~ Compartment)) %>%
    group_by(Region, Date, Compartment, SimID, Source) %>%
    dplyr::summarize(Cases=sum(Cases)) %>%
    ungroup()


plot_order = c('Reported hospital deaths','Reported cumulative hospital deaths', 'Reported census hospital beds occupied',
  'Total incident deaths', 'Reported cumulative deaths', 'new_hospitalizations', 
  'Daily new infections', 'Total infectious', 'Fraction recovered', 
  'HFR','IHR', 'Fraction of deaths outside the hospital',
  'IFR', 'Transmission rate', 'R(t)',
  'Hospital duration (recover)', 'Hospital duration (death)', 'Fraction of non-hospitalized infections that die')

compartment_order = c('EMResource', 'Simulation', 'Simulation (reported)', 'Linelist', 'CLI')

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
    mutate(Compartment='Reported cumulative hospital deaths') %>%
    mutate(Source='Simulation (reported)')

cumulative_all_deaths = plotout_noage %>%  
    filter(Compartment=="Total incident deaths", Source=='Simulation (reported)') %>%
    group_by(SimID) %>%
    arrange(Date) %>%
    mutate(Cases=cumsum(Cases)) %>%
    ungroup() %>%
    mutate(Compartment='Reported cumulative deaths') %>%
    mutate(Source = 'Simulation (reported)')

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
phi = unlogit(pars[[paste0('phi_scale_', region_to_test)]], pars$IFR_constraint,0)
IHR_init = (pars$IFR_constraint - phi) / (HFR_init - phi)
paras = list(sigma = 1/3.5,
             zeta_s = 1/1.5,
             mu_m = 1/15,
             gamma_m = 1/3.5,
             zeta_h = 1/4.2,
             phi = phi)
r0pop = read.csv(covid_get_path('Data/covid_region_populations.csv')) %>% 
    summarize(POPULATION=sum(POPULATION))

source('R0_functions.R')
R0_val = get_R0(region_cons=region_to_test, 
       beta_value=beta_init,
       pop=r0pop,
       paras=paras,
       t=47,
       IHR =IHR_init)

Rt_scaling_constant = R0_val / beta_init
print(paste(R0_val, Rt_scaling_constant, beta_init, phi, HFR_init, IHR_init))
plotout_noage = bind_rows(plotout_noage %>% filter(Compartment != "Fraction recovered"), frac_recovered)
Rt = plotout_noage %>%
  filter(Compartment == 'Transmission rate' | Compartment == 'Fraction recovered') %>%
  spread(Compartment, Cases) %>%
  mutate(Rt = Rt_scaling_constant * `Transmission rate` * (1 - `Fraction recovered`)) %>%
  select(-`Transmission rate`, -`Fraction recovered`) %>%
  rename(Cases=Rt) %>%
  mutate(Compartment='R(t)')
plotout_noage = bind_rows(plotout_noage, Rt)

rt_model = Rt %>%
  select(-Time) %>%
  drop_na() %>%
  group_by(Date) %>%
  summarize(rt_model_med=quantile(Cases, 0.5, type=3),
            rt_model_lower=quantile(Cases, 0.025, type=3),
            rt_model_upper=quantile(Cases, 0.975, type=3)) %>%
  mutate(geography_modeled = sprintf('covidregion_%s', region_to_test),
         scenario_name = 'baseline') %>%
  ungroup() %>%
  rename(date=Date)

write.csv(rt_model, sprintf('%s/model_rt_%s_%s.csv', human_output_dir, region_to_test, Sys.Date()-1), row.names=F)

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
phidf = data.frame(Date = as.Date('2020-01-14') + times, Compartment='Fraction of non-hospitalized infections that die', 
  Cases=unlogit(pars[[paste0('phi_scale_', region_to_test)]], pars[['IFR_constraint']], 0), Source='Simulation')

plotout_noage = bind_rows(plotout_noage, gamma_df) %>% bind_rows(mu_df) %>% bind_rows(phidf)

print("Plotting")
ggplot(plotout_noage %>% filter(Date >= project_zoom), aes(x=Date, y=Cases, color=Source, fill=Source)) +
   stat_summary(fun=function(z){quantile(z,0.5, type=3)}, geom="line", size=0.75) +
   stat_summary(fun.min=function(z){quantile(z,0.025, type=3)}, fun.max=function(z){quantile(z,0.975, type=3)}, geom="ribbon", alpha=0.3, color=NA) +
   geom_point(plotting_data %>% filter(Date >= project_zoom), mapping=aes(x=Date, y=Cases, fill=Source), size=0.9, alpha=0.7, pch=21) +
   facet_wrap(~factor(Compartment, levels=plot_order), scales='free_y', ncol=3) +
   theme_bw() +
   ylab('') +
   scale_x_date(date_breaks = "months", date_labels="%b") +
   ggtitle(sprintf("Region %s\nR0: %s\nInitial prevalence: %s\nphi: %s", 
    region_to_test, R0_val, round(popfinal[[region_to_test]] * pars[[paste0('num_init_', region_to_test)]]), mean(phidf$Cases))) +
   scale_color_brewer(palette='Set1', limits=compartment_order) +
   scale_fill_brewer(palette='Set1', limits=compartment_order) +
   geom_vline(xintercept=Sys.Date() - 1)
ggsave(sprintf('%s/projection_%s_%s.png', human_output_dir, model_name, region_to_test), width=14, height=8)


