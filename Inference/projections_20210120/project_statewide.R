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

## read in plotting data
plotting_data = read.csv(covid_get_path(emr_data_file)) %>%
  filter(covid_region %in% 1:11) %>%
  select(time, ObsHosp_1,ObsICU_1) %>%
  mutate(Date = as.Date('2020-01-14') + time) %>%
  group_by(Date) %>%
  summarize(ObsHosp_1 = sum(ObsHosp_1), ObsICU_1 = sum(ObsICU_1)) %>%
  ungroup() %>%
  arrange(Date) %>%
  select(Date, ObsHosp_1, ObsICU_1) %>%
  gather(Compartment, Cases, ObsHosp_1,  ObsICU_1) %>%
  mutate(Compartment = case_when(
                            (Compartment == 'ObsHosp_1') ~ 'Reported census hospital beds occupied',
                            (Compartment == 'ObsICU_1') ~ 'Census ICU beds occupied'),
    Source = 'EMResource') %>%
  select(Date, Compartment, Cases, Source)

ll_deaths = read.csv(covid_get_path(ll_file)) %>%
  mutate(ll_deaths = ll_deaths - ll_hosp_deaths) %>%
  mutate(date = as.Date(date)) %>%
  rename(Date = date, Cases = ll_deaths) %>%
  mutate(Compartment = 'Reported non-hospital deaths') %>%
  mutate(Source='Linelist') %>%
  group_by(Date, Compartment, Source) %>%
  summarize(Cases=sum(Cases)) %>%
  ungroup()

cumulative_ll = ll_deaths %>%
  arrange(Date) %>%
  mutate(Cases = cumsum(Cases),
    Compartment = 'Reported cumulative non-hosp deaths')

ll_hosp_deaths = read.csv(covid_get_path(ll_file)) %>%
  mutate(date = as.Date(date)) %>%
  rename(Date = date, Cases = ll_hosp_deaths) %>%
  mutate(Compartment = 'Reported hospital deaths') %>%
  mutate(Source='Linelist') %>%
  group_by(Date, Compartment, Source) %>%
  summarize(Cases=sum(Cases)) %>%
  ungroup()
cumulative_hosp_ll = ll_hosp_deaths %>%
  arrange(Date) %>%
  mutate(Cases = cumsum(Cases),
    Compartment = 'Reported cumulative hospital deaths')

cli_data = read.csv(covid_get_path(cli_file)) %>%
  mutate(Date=as.Date(date),
    Compartment='new_hospitalizations',
    Source='CLI') %>%
  rename(Cases = cli) %>%
  select(Date, Cases, Compartment, Source) %>%
  group_by(Date, Compartment, Source) %>%
  summarize(Cases=sum(Cases)) %>%
  ungroup()%>%
  filter(Date < (max(Date) - 3))

plotting_data = plotting_data %>%
  bind_rows(ll_deaths) %>%
  bind_rows(cli_data) %>%
  bind_rows(cumulative_ll) %>%
  bind_rows(ll_hosp_deaths) %>%
  bind_rows(cumulative_hosp_ll)


print("Writing output")
## Output civis format

foreach (r=1:11, .combine='rbind') %do%{
  p = read.csv(sprintf("%s/projections_%s_%s.csv", human_output_dir, model_name, r))
  p
} ->alloutputs

alloutputs = alloutputs %>%
    mutate(Date=as.Date(date)) %>%
    select(-date) %>%
    group_by(SimID, Date, Compartment) %>%
    summarize(Cases = sum(Cases)) %>%
    ungroup() %>%
    mutate(Compartment = as.character(Compartment)) %>%
    mutate(Time=as.numeric(Date - as.Date('2020-01-14'))) %>%
    mutate(Region = 'Illinois')

civis_proj = alloutputs %>% rename(date=Date)


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
    
    cases = civis_proj %>%
        filter(date >= as.Date('2020-03-13')) %>%
        filter(Compartment == comp) %>%
        group_by(date) %>%
        dplyr::summarize(med=quantile(Cases, 0.5, type=3),
                         lower=quantile(Cases, 0.025, type=3),
                         upper=quantile(Cases, 0.975, type=3)) %>%
        mutate(geography_modeled = 'illinois',
               scenario_name = 'baseline') %>%
        ungroup() %>%
        select(date, geography_modeled, scenario_name, med, lower, upper)
    names(cases) = c('date','geography_modeled','scenario_name',paste0(comp_name, c('median','lower','upper')))
    cases
} -> tempframe
finaldf = tempframe %>% reduce(full_join, by = c("date", "scenario_name", "geography_modeled"))
write.csv(finaldf, sprintf('%s/uchicago_covidregion_%s_%s.csv', human_output_dir, 'illinois', Sys.Date()-1), row.names=F)



compartments = c('ObsHospDeaths', 'ObsDeaths', "ObsHosp", "ObsICU", 'Inc', "Total infectious",  "R", 'new_hospitalizations')
new_compartment_name_map = c('Reported hospital deaths', 'Reported non-hospital deaths', 'Reported census hospital beds occupied', 'Census ICU beds occupied',
  'Daily new infections', 'Total infectious', 'Fraction recovered', 'new_hospitalizations')
names(new_compartment_name_map) = compartments
plot_order = c('Reported hospital deaths','Reported cumulative hospital deaths', 'Reported census hospital beds occupied', 'Census ICU beds occupied',
  'Reported non-hospital deaths', 'Reported cumulative non-hosp deaths', 'new_hospitalizations', 
  'Daily new infections', 'Total infectious', 'Fraction recovered')
compartment_order = c('EMResource', 'Simulation', 'Simulation (reported)', 'Linelist', 'CLI')

plotout = alloutputs %>%
  filter(Compartment %in% compartments) %>%
  mutate(Compartment = new_compartment_name_map [Compartment],
    Source = ifelse(Compartment %in% c('ObsHospDeaths', 'ObsDeaths', 'ObsHosp', "ObsICU"), 'Simulation (reported)', 'Simulation'),
    Cases = ifelse(Compartment == 'Fraction recovered', Cases/12821497, Cases))

cumulative_hosp = get_cumulative(plotout, "Reported hospital deaths", 'Reported cumulative hospital deaths', 'Simulation (reported)')
cumulative_all_deaths = get_cumulative(plotout, "Reported non-hospital deaths", 'Reported cumulative non-hosp deaths', 'Simulation (reported)')

plotout = bind_rows(plotout, cumulative_hosp) %>%
  bind_rows(cumulative_all_deaths)

ggplot(plotout, aes(x=Date, y=Cases, color=Source, fill=Source)) +
     stat_summary(fun=function(z){quantile(z,0.5,type=3)}, geom="line", size=0.75) +
     stat_summary(fun.min=function(z){quantile(z,0.025, type=3)}, fun.max=function(z){quantile(z,0.975, type=3)}, geom="ribbon", alpha=0.3, color=NA) +
     geom_point(plotting_data, mapping=aes(x=Date, y=Cases, fill=Source), size=0.9, alpha=0.7, pch=21) +
     facet_wrap(~factor(Compartment, levels=plot_order), scales='free_y', ncol=3) +
     theme_bw() +
     ylab('') +
     scale_x_date(date_breaks = "months", date_labels="%b") +
     scale_color_brewer(palette='Set1', limits=compartment_order) +
     scale_fill_brewer(palette='Set1', limits=compartment_order)
ggsave(sprintf('%s/projection_%s_%s.png', human_output_dir, model_name, 'illinois'), width=14, height=10)
