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
library(zoo)

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
region_to_test=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
stopifnot(region_to_test %in% c(-1, 1:11))

## read in plotting data
plotting_data = read.csv(covid_get_path(emr_data_file)) %>%
  filter(covid_region == region_to_test) %>%
  select(time, ObsHosp_1, ObsHospDeaths_1, ObsICU_1) %>%
  mutate(Date = as.Date('2020-01-14') + time) %>%
  arrange(Date) %>%
  select(Date, ObsHosp_1, ObsICU_1) %>%
  gather(Compartment, Cases, ObsHosp_1, ObsICU_1) %>%
  mutate(Compartment = case_when(
                            (Compartment == 'ObsHosp_1') ~ 'Reported census hospital beds occupied',
                            (Compartment == 'ObsICU_1') ~ 'Census ICU beds occupied'),
    Source = 'EMResource') %>%
  select(Date, Compartment, Cases, Source) %>%
  mutate(Date = Date)

ll_deaths = read.csv(covid_get_path(ll_file)) %>%
  mutate(ll_deaths = ll_deaths - ll_hosp_deaths) %>%
  filter(covid_region == region_to_test) %>%
  mutate(date = as.Date(date)) %>%
  rename(Date = date, Cases = ll_deaths) %>%
  mutate(Compartment = 'Reported non-hospital deaths') %>%
  mutate(Source='Linelist')

cumulative_ll = ll_deaths %>%
  arrange(Date) %>%
  mutate(Cases = cumsum(Cases),
    Compartment = 'Reported cumulative non-hosp deaths')

ll_hosp_deaths = read.csv(covid_get_path(ll_file)) %>%
  filter(covid_region == region_to_test) %>%
  mutate(date = as.Date(date)) %>%
  rename(Date = date, Cases = ll_hosp_deaths) %>%
  mutate(Compartment = 'Reported hospital deaths') %>%
  mutate(Source='Linelist')
cumulative_hosp_ll = ll_hosp_deaths %>%
  arrange(Date) %>%
  mutate(Cases = cumsum(Cases),
    Compartment = 'Reported cumulative hospital deaths')


cli_data = read.csv(covid_get_path(cli_file)) %>%
  mutate(Date=as.Date(date),
    Compartment='new_hospitalizations',
    Source='CLI') %>%
  rename(Cases = cli) %>%
  filter(covid_region == region_to_test) %>%
  select(Date, Cases, Compartment, Source) %>%
  filter(Date < (max(Date) - 3))

plotting_data = plotting_data %>%
  bind_rows(ll_deaths) %>%
  bind_rows(cli_data) %>%
  bind_rows(cumulative_ll) %>%
  bind_rows(ll_hosp_deaths) %>%
  bind_rows(cumulative_hosp_ll)

pars$region_to_test = region_to_test
regtemp = region_to_test
pars$tmin = simstart
pars$tmax = project_end

pops = population_list %>%
    group_by(covid_region) %>%
    summarize(pop=sum(POPULATION)) %>%
    ungroup() 
popfinal = pops %>% select(pop) %>% unlist()
names(popfinal) = pops$covid_region


  sim_0 = 0
  slice_df = read.csv(sprintf("%s/num_sims_%s_%s.csv", output_dir, model_name, region_to_test)) %>%
    filter(replicates > 0) %>%
    select(file, replicates, loglik)

  AIC = max(slice_df[['loglik']])

  foreach (r=1:nrow(slice_df), .combine='rbind') %do%{

    pf = readRDS(slice_df[r, 'file'])
    n_sim_project = slice_df[r, 'replicates']
    new_pars = pf@params

    ## Simulate
    sim_temp = simulate_pomp_covid(
      n_regions = pars$n_regions,
      n_age_groups = n_age_groups,
      nsim=n_sim_project, 
      input_params = pars,
      delta_t = deltaT,
      population_list=population_list,
      death_reporting_covar=fraction_underreported,
      emr_report_covar=emr_report_covar,
      rprocess_Csnippet = rprocess_snippet,
      rinit_Csnippet = rinit_snippet,
      global_Csnippet = global_snippet,
      rmeasure_Csnippet = rmeasure_snippet,
      obsnames = c(paste0('ObsHosp_', 1), paste0('ObsHospDeaths_',1), paste0('ObsDeaths_', 1), paste0('ObsICU_', 1)),
      use.mle.params = T,
      mle.params=new_pars
      )
    sim_temp$raw_simulation_output = sim_temp$raw_simulation_output %>%
      mutate(.id = as.numeric(.id), .id = .id + sim_0)
    sim_0 = sim_0 + n_sim_project
    sim_temp$raw_simulation_output
  } -> sim_full

  ## Process output

alloutputs = process_pomp_covid_output(sim_full, agg_regions=F)


# Make basic plots
compartments_to_plot = c('ObsHospDeaths', 'ObsDeaths', "ObsHosp", "ObsICU", 'Inc', "Total infectious",  "R", 'new_hospitalizations', 
      'transmissionRate', 'HFRtrack','IHRtrack')
new_compartment_name_map = c('Reported hospital deaths', 'Reported non-hospital deaths', 'Reported census hospital beds occupied', 'Census ICU beds occupied',
  'Daily new infections', 'Total infectious', 'Fraction recovered', 'new_hospitalizations', 'Transmission rate',
  'HFR','IHR')
names(new_compartment_name_map) = compartments_to_plot
plot_order = c('Reported hospital deaths','Reported cumulative hospital deaths', 'Reported census hospital beds occupied', 'Census ICU beds occupied',
  'Reported non-hospital deaths', 'Reported cumulative non-hosp deaths', 'new_hospitalizations', 
  'Daily new infections', 'Total infectious', 'Fraction recovered', 
  'HFR','IHR',
  'IFR', 'Transmission rate', 'R(t)',
  'Hospital duration (recover)', 'Hospital duration (death)', 'ICU duration', 'ICU fraction')
compartment_order = c('EMResource', 'Simulation', 'Simulation (reported)', 'Linelist', 'CLI')
plot_generic(alloutputs, 
  plotting_data, 
  new_pars, 
  compartments_to_plot, 
  new_compartment_name_map, 
  plot_order, 
  compartment_order,
  AIC, 
  popfinal)


# Make projection outputs
projections = alloutputs$plotting_output %>%  
    mutate(date=as.Date('2020-01-14') + Time) %>%
    select(-Time)
projections_nonicu = projections %>%
  filter(Compartment %in% c('IH', 'IC')) %>% 
  spread(Compartment, Cases) %>%
  mutate(Cases = IH - IC,
    Cases=ifelse(Cases<0, 0, Cases)) %>%
  mutate(Compartment='IH_nonicu') %>%
  ungroup() %>%
  select(-IH, -IC)
projections = bind_rows(projections, projections_nonicu)
proj_deaths = projections %>%
  filter(Compartment %in% c('ObsDeaths', 'ObsHospDeaths')) %>% 
  spread(Compartment, Cases) %>%
  mutate(Cases = ObsDeaths + ObsHospDeaths) %>%
  mutate(Compartment='ObsTotDeaths') %>%
  ungroup() %>%
  select(-ObsDeaths, -ObsHospDeaths)
projections = projections %>%
  bind_rows(proj_deaths)
write.csv(projections, sprintf("%s/projections_%s_%s.csv", human_output_dir, model_name, region_to_test), row.names=F)

compartments = c(`Total infectious`='cases_',
  Inc='cases_new_',
  new_deaths='deaths_',
  ObsTotDeaths='deaths_det_',
  IH_nonicu='hosp_bed_',
  R='recovered_',
  IC='icu_',
  IH='total_hosp_')
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
        filter(date <= d & date>=Sys.Date()-1) %>%
        filter(Compartment == 'IH') %>%
        group_by(SimID) %>%
        dplyr::summarize(exceeded = any(Cases >= beds)) %>%
        select(exceeded) %>%
        unlist(use.names=F)
    hb_template[r, 'percent_of_simulations_that_exceed'] = sum(temp)/nsim * 100
}

ic_template = read.csv(covid_get_path(hosp_capacity_file)) %>%
  filter(geography_modeled == paste0('covidregion_', region_to_test)) %>%
  filter(resource_type=='icu_availforcovid')
nsim = length(unique(projections$SimID))
for (r in 1:nrow(ic_template)){
    beds = as.numeric(ic_template[r, 'avg_resource_available'])
    d = as.Date(ic_template[r, 'date_window_upper_bound'])

    temp = projections %>% 
        filter(date <= d & date>=(Sys.Date()-1)) %>%
        filter(Compartment == 'IC') %>%
        group_by(SimID) %>%
        dplyr::summarize(exceeded = any(Cases >= beds)) %>%
        select(exceeded) %>%
        unlist(use.names=F)

    ic_template[r, 'percent_of_simulations_that_exceed'] = sum(temp)/nsim * 100
}

hb_template = bind_rows(ic_template, hb_template)

write.csv(hb_template %>% mutate(scenario_name='baseline'), sprintf('%s/uchicago_hospitaloverflow_covidregion_%s_%s.csv', human_output_dir, region_to_test, Sys.Date()-1), row.names=F)


