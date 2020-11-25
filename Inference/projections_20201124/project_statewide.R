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
  select(time, ObsHosp_1, ObsHospDeaths_1) %>%
  mutate(Date = as.Date('2020-01-14') + time) %>%
  group_by(Date) %>%
  summarize(ObsHosp_1 = sum(ObsHosp_1), ObsHospDeaths_1 = sum(ObsHospDeaths_1)) %>%
  ungroup() %>%
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
  mutate(date = as.Date(date)) %>%
  rename(Date = date, Cases = ll_deaths) %>%
  mutate(Compartment = 'Total incident deaths') %>%
  mutate(Source='Linelist') %>%
  group_by(Date, Compartment, Source) %>%
  summarize(Cases=sum(Cases)) %>%
  ungroup()

cumulative_ll = ll_deaths %>%
  arrange(Date) %>%
  mutate(Cases = cumsum(Cases),
    Compartment = 'Reported cumulative deaths')

cli_data = read.csv(covid_get_path(cli_file)) %>%
  mutate(Date=as.Date(date),
    Compartment='new_hospitalizations',
    Source='CLI') %>%
  rename(Cases = cli) %>%
  select(Date, Cases, Compartment, Source) %>%
  group_by(Date, Compartment, Source) %>%
  summarize(Cases=sum(Cases)) %>%
  ungroup()

plotting_data = plotting_data %>%
  bind_rows(ll_deaths) %>%
  bind_rows(cli_data) %>%
  bind_rows(cumulative_ll)


print("Writing output")
## Output civis format

alloutputs = read.csv(sprintf("%s/projections_%s_%s.csv", human_output_dir, model_name, 'illinois')) %>%
    mutate(Date=as.Date(date)) %>%
    select(-date, -X) %>%
    mutate(Compartment = as.character(Compartment)) %>%
    mutate(Time=as.numeric(Date - as.Date('2020-01-14')))

print(head(alloutputs))
plotout_noage = alloutputs %>% ungroup() %>%    
    filter(Compartment %in% c('ObsHospDeaths', 'ObsDeaths', "ObsHosp", 'Inc', "Total infectious",  "R", 'new_hospitalizations', 'new_deaths')) %>%
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
                            T ~ Compartment)) %>%
    group_by(Date, Compartment, SimID, Source) %>%
    dplyr::summarize(Cases=sum(Cases)) %>%
    ungroup()


plot_order = c('Reported hospital deaths','Reported cumulative hospital deaths', 'Reported census hospital beds occupied',
  'Total incident deaths', 'Reported cumulative deaths', 'new_hospitalizations', 
  'Daily new infections', 'Total infectious', 'Fraction recovered', 'IFR')

compartment_order = c('EMResource', 'Simulation', 'Simulation (reported)', 'Linelist', 'CLI')


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

D= alloutputs %>% ungroup() %>%  
    filter(Compartment=="D") %>%
    group_by(SimID, Time) %>%
    summarize(Dead=sum(Cases)) %>%
    ungroup() %>%
    mutate(Time = Time - 19)

R= alloutputs %>% ungroup() %>%  
    filter(Compartment=="R") %>%
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
  mutate(Cases = Cases / 12821497)

plotout_noage = bind_rows(plotout_noage %>% filter(Compartment != "Fraction recovered"), frac_recovered)

print("Plotting")
ggplot(plotout_noage %>% filter(Date >= project_zoom), aes(x=Date, y=Cases, color=Source, fill=Source)) +
   stat_summary(fun=function(z){quantile(z,0.5, type=3), geom="line", size=0.75) +
   stat_summary(fun.min=function(z){quantile(z,0.025, type=3)}, fun.max=function(z){quantile(z,0.975, type=3)}, geom="ribbon", alpha=0.3, color=NA) +
   geom_point(plotting_data %>% filter(Date >= project_zoom), mapping=aes(x=Date, y=Cases, fill=Source), size=0.9, alpha=0.7, pch=21) +
   facet_wrap(~factor(Compartment, levels=plot_order), scales='free_y', ncol=3) +
   theme_bw() +
   ylab('') +
   scale_x_date(date_breaks = "months", date_labels="%b") +
   ggtitle("Illinois") +
   scale_color_brewer(palette='Set1', limits=compartment_order) +
   scale_fill_brewer(palette='Set1', limits=compartment_order) +
   geom_vline(xintercept=Sys.Date() - 1)
ggsave(sprintf('%s/projection_%s_%s.png', human_output_dir, model_name, "Illinois"), width=14, height=8)


