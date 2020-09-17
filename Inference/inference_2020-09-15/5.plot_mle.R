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

source('./input_file_specification.R')
default_par_file = './final_mle_pars.csv'

# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1

# Reading input args
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
region_to_test=-1

# Read in inference functions
source(covid_get_path(inference_file))
source('set_up_covariates_and_data.R')
pars$region_to_test = region_to_test
pars$tmin = simstart

## Reformat civis data
region_order = as.character(1:pars$n_regions)
region_plot_order=c('COVID 1','COVID 2', 'COVID 3 and 6', 'COVID 7-10', 'COVID 4-5', 'COVID 11', 'Illinois')

df = df %>% mutate(restore_region = region_plot_order[restore_region])

IL_civis = df %>%
    group_by(date, source) %>%
    summarize(emr_deaths = sum(emr_deaths),
              nonhosp_deaths = sum(nonhosp_deaths),
              confirmed_covid_icu = sum(confirmed_covid_icu),
              total_deaths=sum(total_deaths),
              incident_hospitalizations=sum(new_ll_hospitalizations),
              confirmed_covid_nonicu=sum(covid_non_icu)) %>%
        ungroup() %>%
    mutate(restore_region = 'Illinois')
    
civis_data = bind_rows(df, IL_civis) %>%  
  mutate(Date = date) %>%
  select(-date)
print(pars)
## Simulate
sim_full = simulate_pomp_covid(
  n_regions = pars$n_regions,
  n_age_groups = n_age_groups,
  nsim=50, 
  input_params = pars,
  delta_t = deltaT,
  population=population_list,
  beta_scales=beta_scales,
  beta_covar=beta_covariate_column,
  contacts=pomp_contacts,
  frac_underreported=fraction_underreported,
  rprocess_Csnippet = rprocess_snippet,
  rinit_Csnippet = rinit_snippet,
  rmeasure_Csnippet=rmeasure_snippet,
  obsnames =observed_names
)
print(pars)
## Process output
alloutputs = process_pomp_covid_output(sim_full, agg_regions=F)
write.csv(sim_full$raw_simulation_output, 'raw_output.csv',  row.names=F)
pop_totals = population_list %>% group_by(region) %>%
  summarize(POPULATION = sum(POPULATION)) %>%
  select(POPULATION) %>%
  unlist(use.names=F)
pop_totals = c(pop_totals, sum(pop_totals))
names(pop_totals) = region_plot_order

plotout = alloutputs$plotting_output %>% ungroup() %>%    
    mutate(Date=as.Date('2020-01-14') + Time,
      restore_region = region_plot_order[as.numeric(Region)]) 
df_infections_summary <- plotout %>%
    ungroup() %>% 
    group_by(Time, Compartment, Region) %>%
    filter(Compartment %in% c('E', 'A', 'P', 'IS', 'IM', 'IH', 'IC')) %>%
    group_by(SimID, Date, Region) %>%
    arrange(Time) %>%
    summarise(Compartment = 'Prevalence',
              Cases = sum(Cases)) %>%
  ungroup()

statewide = plotout %>%
    group_by(Date, SimID, Compartment) %>%
    summarize(Cases=sum(Cases)) %>%
    mutate(Region='0', restore_region='Illinois') %>%
    ungroup()
plotout = bind_rows(plotout, statewide)

## Plot
write.csv(plotout, './mle_plotting_data.csv', row.names=F)



foreach(rnum=1:7) %do%{
    region_to_plot = region_plot_order[rnum]

    civisplot = civis_data %>%
        filter(restore_region==region_to_plot) %>% 
        mutate(deaths=total_deaths) %>%
        select(-covid_non_icu, -total_deaths, -new_ll_hospitalizations) %>%
        gather(Compartment, 
          Cases, 
          deaths, 
          confirmed_covid_icu, 
          confirmed_covid_nonicu, 
          emr_deaths, 
          nonhosp_deaths, 
          incident_hospitalizations) %>% 
        mutate(source=case_when((Compartment=='confirmed_covid_icu') ~ 'emresource',
                                (Compartment=='deaths') ~ source,
                                (Compartment=='emr_deaths') ~ 'emresource',
                                (Compartment=='nonhosp_deaths') ~ 'IDPH line list',
                                (Compartment=='incident_hospitalizations') ~ 'IDPH line list',
                                (Compartment=='confirmed_covid_nonicu') ~ 'emresource')) %>%
        mutate(Compartment=case_when((Compartment=='deaths')~'All reported deaths',
                                     (Compartment=='confirmed_covid_icu')~'Census ICU beds occupied',
                                     (Compartment=='nonhosp_deaths')~'Non-hospital deaths',
                                     (Compartment=='emr_deaths')~'Hospitalized deaths',
                                     (Compartment=='incident_hospitalizations') ~ 'Incident hospitalizations',
                                     (Compartment=='confirmed_covid_nonicu') ~ 'Census non-ICU beds occupied' ))
    print(tail(civisplot %>% filter(Compartment =='Census non-ICU beds occupied')))
    plotout %>% 
        filter(Compartment %in% c('Reported deaths', 'ObsICU', 'nNHD', 'new_hospitalizations', 'ObsHospDeaths', 'ObsHosp')) %>% 
        filter(restore_region==region_to_plot) %>%
        mutate(Compartment=case_when((Compartment=='Reported deaths')~'All reported deaths',
                                     (Compartment=='ObsICU')~'Census ICU beds occupied',
                                     (Compartment=='nNHD')~'Non-hospital deaths',
                                     (Compartment=='ObsHospDeaths')~'Hospitalized deaths',
                                     (Compartment=='new_hospitalizations') ~ 'Incident hospitalizations',
                                     (Compartment=='ObsHosp') ~ 'Census non-ICU beds occupied')) %>% 
        spread(Compartment, Cases) %>% 
        mutate(`Non-hospital deaths` = `All reported deaths` - `Hospitalized deaths`,
              `Non-hospital deaths` = if_else(`Non-hospital deaths` < 0, 0, `Non-hospital deaths`)) %>%
        gather(Compartment, Cases, 
          'All reported deaths', 'Census ICU beds occupied', 'Non-hospital deaths', 'Hospitalized deaths', 'Incident hospitalizations', 'Census non-ICU beds occupied') %>%
        filter(Compartment != 'Incident hospitalizations') %>% 

        ggplot(aes(x=Date, y=Cases)) + 
        stat_summary(fun=median, geom="line", color="black", size=1) +
        stat_summary(fun.min=function(z){quantile(z,0.025)}, fun.max=function(z){quantile(z,0.975)}, geom="ribbon", color="black", alpha=0.1) +
        geom_point(civisplot %>% filter(Compartment != 'Incident hospitalizations'), mapping=aes(x=Date, y=Cases, fill=source), color='black', pch=21, size=2, alpha=0.5) +
        facet_wrap(~Compartment, scales='free_y', ncol=5) +
        ylab(region_to_plot) + 
        theme_bw() +
        scale_fill_brewer(palette='Dark2') -> p
    p
} -> plots

g = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], nrow=7)
ggsave('mle_fit.png', g, width=20, height=16)




order_to_plot_states = c('Prevalence',
                         'Incidence',
                         'Recovered',
                         'Susceptible',
                         'Incident deaths',
                         'Non-ICU hospital beds occupied',
                         'ICU beds occupied')

foreach(rnum=1:7) %do%{
    region_to_plot = region_plot_order[rnum]
    popregion = pop_totals[region_to_plot]
plotout %>% 
    filter(Compartment %in% c('Incidence', 'Prevalence','IH','IC', 'nD', 'R', 'S')) %>% 
    mutate(Compartment = case_when((Compartment=='Incidence') ~ 'Incidence',
                                   (Compartment=='Prevalence') ~ 'Prevalence',
                                   (Compartment == 'IH') ~ 'Non-ICU hospital beds occupied',
                                   (Compartment == 'IC') ~ 'ICU beds occupied',
                                   (Compartment == 'nD') ~ 'Incident deaths',
                                   (Compartment == 'R') ~ 'Recovered',
                                   (Compartment == 'S') ~ 'Susceptible')) %>%
    filter(restore_region==region_to_plot) %>%
    mutate(Cases = case_when((Compartment %in% c('Prevalence','Recovered', 'Susceptible')) ~ Cases / popregion,
                            (!Compartment %in% c('Prevalence','Recovered', 'Susceptible')) ~ Cases)) %>%
    mutate(Compartment = factor(Compartment, levels=order_to_plot_states))  -> incidence_plots

incidence_plots %>%
    filter(Compartment!='Susceptible') %>%
    ggplot(aes(x=Date, y=Cases)) + 
    stat_summary(fun.ymin=function(z){quantile(z,0.025)}, fun.ymax=function(z){quantile(z,0.975)}, geom="ribbon", color="black", alpha=0.1) +
    stat_summary(fun.y=median, geom="line", color="black", size=1) +
    facet_wrap(~Compartment, scales='free_y', ncol=7) +
    ylab(region_to_plot) + 
    theme_bw()->p
p
} -> plots

g = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], nrow=7)
ggsave('mle_latent_states.png', g, width=32, height=20)