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

region_order = c('northcentral','central','northeast','southern', 'chicago')
civis_data = read.csv(data_filename)  %>% 
    mutate(Date=as.Date(date), total_deaths = hosp_deaths+nonhosp_deaths) %>%
    filter(Date >= as.Date('2020-03-16'))

LL_end = as.Date('2020-06-13')

get_idph = function(date, region){
    idph[which(idph$Date==date & idph$restore_region==region), 'total_deaths']
}

# Add smoothed total deaths to total deaths
idph = read.csv(idph_filename) %>% 
    mutate(Date=as.Date(test_date), restore_region=region, total_deaths=inc_deaths) %>%
    select(Date, restore_region, total_deaths)

civis_data = civis_data %>% mutate(
    source=case_when((is.na(total_deaths)) ~ 'Public',
                                (!is.na(total_deaths))~ 'IDPH line list'),
    total_deaths=case_when((is.na(total_deaths)) ~ as.numeric(mapply(get_idph, Date, restore_region)),
                                (!is.na(total_deaths))~ as.numeric(total_deaths)),
    incident_hospitalizations=new_ll_hospitalizations,
    confirmed_covid_nonicu=covid_non_icu
                                )
IL_civis = civis_data %>%
    filter(restore_region != 'chicago') %>%
    group_by(date, Date, source, source_hosp_deaths) %>%
    summarize(emr_deaths = sum(emr_deaths),
              nonhosp_deaths = sum(nonhosp_deaths),
              confirmed_covid_icu = sum(confirmed_covid_icu),
              total_deaths=sum(total_deaths),
              incident_hospitalizations=sum(new_ll_hospitalizations),
              confirmed_covid_nonicu=sum(covid_non_icu)) %>%
        ungroup() %>%
    mutate(restore_region = 'Illinois')
    
civis_data = bind_rows(civis_data, IL_civis) %>%
    mutate(confirmed_covid_nonicu = if_else(Date < as.Date('2020-05-06'), NA_integer_, confirmed_covid_nonicu))

## Simulate
sim_full = simulate_pomp_covid(
  n_regions = pars$n_regions,
  n_age_groups = n_age_groups,
  nsim=5, 
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

## Process output
alloutputs = process_pomp_covid_output(sim_full, agg_regions=F)

pop_totals = as.numeric(lapply(FUN=sum, population_list))
chicago_pop = pop_totals[5]
pop_totals = pop_totals[1:4]
pop_totals = c(pop_totals, chicago_pop, sum(pop_totals))
names(pop_totals) = c(region_order[1:4], 'chicago', 'Illinois')

plotout = alloutputs$plotting_output %>% ungroup() %>%    
    mutate(Date=as.Date('2020-01-14') + Time)
df_infections_summary <- plotout %>%
    ungroup() %>% 
    group_by(Time, Compartment, Region) %>%
    filter(Compartment %in% c('E', 'A', 'P', 'IS', 'IM', 'IH', 'IC')) %>%
    group_by(SimID, Date, Region) %>%
    arrange(Time) %>%
    summarise(Compartment = 'Prevalence',
              Cases = sum(Cases)) %>%
  ungroup()
chicago_out = bind_rows(plotout, df_infections_summary) %>% 
    filter(Region==5) %>%
    mutate(restore_region=region_order[as.numeric(Region)]) %>%
    select(-Time)
plotout = bind_rows(plotout, df_infections_summary) %>% 
    filter(Region!=5) %>% # Remove Chicago
    mutate(restore_region=region_order[as.numeric(Region)]) %>%
    select(-Time)
statewide = plotout %>%
    group_by(Date, SimID, Compartment) %>%
    summarize(Cases=sum(Cases)) %>%
    mutate(Region='0', restore_region='Illinois') %>%
    ungroup()
plotout = bind_rows(plotout, statewide, chicago_out)

## Plot
write.csv(plotout, './_plotting_data.csv', row.names=F)

region_plot_order=c('northeast','southern','central','northcentral', 'chicago', 'Illinois')

foreach(rnum=1:6) %do%{
    region_to_plot = region_plot_order[rnum]

    civisplot = civis_data %>% 
        mutate(deaths=total_deaths) %>%
        gather(Compartment, Cases, deaths, confirmed_covid_icu, confirmed_covid_nonicu, emr_deaths, nonhosp_deaths, incident_hospitalizations) %>% 
        filter(restore_region==region_to_plot) %>%
        mutate(source=case_when((Compartment=='confirmed_covid_icu') ~ 'emresource',
                                (Compartment=='deaths') ~ source,
                                (Compartment=='emr_deaths') ~ 'emresource',
                                (Compartment=='nonhosp_deaths') ~ 'IDPH line list',
                                (Compartment=='incident_hospitalizations') ~ 'IDPH line list',
                                (Compartment=='confirmed_covid_nonicu') ~ 'emresource')) %>%
        mutate(Compartment=case_when((Compartment=='deaths')~'All reported deaths',
                                     (Compartment=='confirmed_covid_icu')~'Confirmed ICU cases',
                                     (Compartment=='nonhosp_deaths')~'Non-hospital deaths',
                                     (Compartment=='emr_deaths')~'Hospitalized deaths',
                                     (Compartment=='incident_hospitalizations') ~ 'Incident hospitalizations',
                                     (Compartment=='confirmed_covid_nonicu') ~ 'Confirmed non-ICU hospital cases' ))

    plotout %>% 
        filter(Compartment %in% c('Reported deaths', 'ObsICU', 'nNHD', 'new_hospitalizations', 'ObsHospDeaths', 'ObsHosp')) %>% 
        filter(restore_region==region_to_plot) %>%
        mutate(Compartment=case_when((Compartment=='Reported deaths')~'All reported deaths',
                                     (Compartment=='ObsICU')~'Confirmed ICU cases',
                                     (Compartment=='nNHD')~'Non-hospital deaths',
                                     (Compartment=='ObsHospDeaths')~'Hospitalized deaths',
                                     (Compartment=='new_hospitalizations') ~ 'Incident hospitalizations',
                                     (Compartment=='ObsHosp') ~ 'Confirmed non-ICU hospital cases')) %>%
        ggplot(aes(x=Date, y=Cases)) + 
        stat_summary(fun=median, geom="line", color="black", size=1) +
        stat_summary(fun.min=function(z){quantile(z,0.025)}, fun.max=function(z){quantile(z,0.975)}, geom="ribbon", color="black", alpha=0.1) +
        geom_point(civisplot, mapping=aes(x=Date, y=Cases, fill=source), color='black', pch=21, size=3, alpha=0.5) +
        geom_smooth(civisplot, mapping=aes(x=Date, y=Cases), color='firebrick', se=F, span=0.75) +
        facet_wrap(~Compartment, scales='free_y', ncol=6) +
        ylab(region_to_plot) + 
        theme_bw() +
        scale_fill_brewer(palette='Dark2') -> p
    p
} -> plots

g = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=6)
ggsave('mle_fit.png', g, width=16, height=16)




region_plot_order=c('northeast','southern','central','northcentral', 'chicago', 'Illinois')
order_to_plot_states = c('Prevalence',
                         'Incidence',
                         'Recovered',
                         'Susceptible',
                         'Incident deaths',
                         'Non-ICU hospital beds occupied',
                         'ICU beds occupied')

foreach(rnum=1:6) %do%{
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

g = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=6)
ggsave('mle_latent_states.png', g, width=32, height=20)