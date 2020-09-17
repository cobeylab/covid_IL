library(dplyr)
library(tidyr)
library(pomp)
library(rmutil)
library(foreach)
library(gridExtra)
library(digest)
library(ggplot2)
library(stringr)
library(data.table)


select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

final_date = as.Date('2020-09-20')
root <- '../../'
source(file.path(root, '_covid_root.R'))

covid_set_root(root)
default_par_file = './default_parameter_values.csv'

source("./input_file_specification_11regions.R") # script to read in necessary files
output_dir = projection_dir
source(covid_get_path(inference_file))
source("./set_up_covariates_and_data_11regions.R")
source(covid_get_path(projection_functions))


# Projection specifications: all projections
reference_date <- as.Date(t_ref) # All numeric times in reference to this date
end_projection_date <- Sys.Date()
start_projection_date <- simstart + reference_date


deltaT = 0.1 # timestep for projections  

regional_aggregation = F # If true, do statewide estimates
initialize=T


civis_data = read.csv(covid_region_data_filename)  %>% 
    mutate(Date=as.Date(date), total_deaths = hosp_deaths+nonhosp_deaths) %>%
    filter(Date >= as.Date('2020-03-16'))

civis_data = civis_data %>% mutate(
    source=case_when((is.na(total_deaths)) ~ 'Public',
                      (!is.na(total_deaths))~ 'IDPH line list'),
    incident_hospitalizations=new_ll_hospitalizations,
    confirmed_covid_nonicu=covid_non_icu
                                )
IL_civis = civis_data %>%
    group_by(date, Date, source, source_hosp_deaths) %>%
    summarize(emr_deaths = sum(emr_deaths),
              nonhosp_deaths = sum(nonhosp_deaths),
              confirmed_covid_icu = sum(confirmed_covid_icu),
              total_deaths=sum(total_deaths),
              incident_hospitalizations=sum(new_ll_hospitalizations),
              confirmed_covid_nonicu=sum(covid_non_icu)) %>%
        ungroup() %>%
    mutate(covid_region = 0)
    
civis_data = bind_rows(civis_data, IL_civis) %>%
    mutate(confirmed_covid_nonicu = if_else(Date < as.Date('2020-05-06'), NA_integer_, confirmed_covid_nonicu),
          confirmed_covid_icu = if_else(Date < as.Date('2020-04-07'), NA_integer_, confirmed_covid_icu))

plotout <- read.csv(sprintf('%s/baseline_all.csv', output_dir)) %>% 
  mutate(Date=as.Date(Date), covid_region = Region) %>%
  select(-X)

### Plotting projections ###
outfile= paste0(output_dir, model_name, '_all.fitplot.png')
outfile2 = paste0(output_dir, model_name, '_all.states.png')

region_plot_order=as.character(0:pars$n_regions) 

foreach(region_to_plot=0:pars$n_regions) %do%{
    civisplot = civis_data %>% 
        mutate(deaths=total_deaths) %>%
        gather(Compartment, Cases, deaths, confirmed_covid_icu, confirmed_covid_nonicu, emr_deaths, nonhosp_deaths, incident_hospitalizations) %>% 
        filter(covid_region==region_to_plot,
               Date < end_projection_date) %>%
        mutate(source=case_when((Compartment=='confirmed_covid_icu') ~ 'emresource',
                                (Compartment=='deaths') ~ source,
                                (Compartment=='emr_deaths') ~ 'emresource',
                                (Compartment=='nonhosp_deaths') ~ 'IDPH line list',
                                (Compartment=='incident_hospitalizations') ~ 'IDPH line list',
                                (Compartment=='confirmed_covid_nonicu') ~ 'emresource')) %>%
        mutate(Compartment=case_when((Compartment=='deaths')~'All reported deaths',
                                     (Compartment=='confirmed_covid_icu')~'Confirmed ICU census',
                                     (Compartment=='nonhosp_deaths')~'Non-hospital deaths',
                                     (Compartment=='emr_deaths')~'Hospitalized deaths',
                                     (Compartment=='incident_hospitalizations') ~ 'Incident hospitalizations',
                                     (Compartment=='confirmed_covid_nonicu') ~ 'Confirmed non-ICU census' )) %>%
        filter(Compartment != "Incident hospitalizations")

    plotout %>% 
        filter(Compartment %in% c('Reported deaths', 'ObsICU', 'nNHD', 'new_hospitalizations', 'ObsHospDeaths', 'ObsHosp')) %>% 
        filter(covid_region==region_to_plot) %>%
        mutate(Compartment=case_when((Compartment=='Reported deaths')~'All reported deaths',
                                     (Compartment=='ObsICU')~'Confirmed ICU census',
                                     (Compartment=='nNHD')~'Non-hospital deaths',
                                     (Compartment=='ObsHospDeaths')~'Hospitalized deaths',
                                     (Compartment=='new_hospitalizations') ~ 'Incident hospitalizations',
                                     (Compartment=='ObsHosp') ~ 'Confirmed non-ICU census')) %>% 
        spread(Compartment, Cases) %>% 
        mutate(`Non-hospital deaths` = `All reported deaths` - `Hospitalized deaths`,
              `Non-hospital deaths` = if_else(`Non-hospital deaths` < 0, 0, as.numeric(`Non-hospital deaths`))) %>%
        gather(Compartment, Cases, 
          'All reported deaths', 'Confirmed ICU census', 'Non-hospital deaths', 'Hospitalized deaths', 'Incident hospitalizations', 'Confirmed non-ICU census') %>%
        drop_na() %>%
        filter(Date <= as.Date(final_date),
          Compartment != "Incident hospitalizations") %>% 

        ggplot(aes(x=Date, y=Cases)) + 
        stat_summary(fun=median, geom="line", color="black", size=1) +
        stat_summary(fun.min=function(z){quantile(z,0.025)}, fun.max=function(z){quantile(z,0.975)}, geom="ribbon", color="black", alpha=0.1) +
        geom_point(civisplot, mapping=aes(x=Date, y=Cases, fill=source), color='black', pch=21, size=2, alpha=0.5) +
        #geom_smooth(civisplot, mapping=aes(x=Date, y=Cases), color='firebrick', se=F, span=0.75) +
        facet_wrap(~Compartment, scales='free_y', ncol=5) +
        ylab(region_to_plot) + 
        theme_bw() +
        scale_fill_brewer(palette='Dark2') -> p

      if (region_to_plot == 0){
        p = p + ylab('Illinois')
      }
    p
} -> plots

g = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], nrow=12)
ggsave(outfile, g, width=25, height=32)



order_to_plot_states = c('Prevalence',
                         'Incidence',
                         'Recovered',
                         'Susceptible',
                         'Rt',
                         'Incident deaths',
                         'Non-ICU hospital beds occupied',
                         'ICU beds occupied')

foreach(region_to_plot=0:pars$n_regions) %do%{
    if (region_to_plot > 0){
          popregion = population_list %>% 
          filter(covid_region == region_to_plot) %>% 
          summarize(total=sum(POPULATION)) %>% 
          unlist(use.names=F)
    } else{
        popregion = sum(population_list$POPULATION)
    }

plotout %>% 
    select(-parset, -SimID) %>%
    filter(Compartment %in% c('Incidence', 'Prevalence','IH','IC', 'nD', 'R', 'S', 'Rt')) %>% 
    mutate(Compartment = case_when((Compartment=='Incidence') ~ 'Incidence',
                                   (Compartment=='Prevalence') ~ 'Prevalence',
                                   (Compartment == 'IH') ~ 'Non-ICU hospital beds occupied',
                                   (Compartment == 'IC') ~ 'ICU beds occupied',
                                   (Compartment == 'nD') ~ 'Incident deaths',
                                   (Compartment == 'R') ~ 'Recovered',
                                   (Compartment == 'S') ~ 'Susceptible',
                                   (Compartment == 'Rt') ~ 'Rt')) %>%
    filter(covid_region==region_to_plot) %>%
    mutate(Cases = case_when((Compartment %in% c('Prevalence','Recovered', 'Susceptible')) ~ as.numeric(Cases / popregion),
                            (!Compartment %in% c('Prevalence','Recovered', 'Susceptible')) ~ as.numeric(Cases))) %>%
    mutate(Compartment = factor(Compartment, levels=order_to_plot_states))  -> incidence_plots

incidence_plots %>%
    filter(Compartment!='Susceptible') %>%
    filter(Date <= as.Date(final_date)) %>% 
    ggplot(aes(x=Date, y=Cases)) + 
    stat_summary(fun.ymin=function(z){quantile(z,0.025)}, fun.ymax=function(z){quantile(z,0.975)}, geom="ribbon", color="black", alpha=0.1) +
    stat_summary(fun.y=median, geom="line", color="black", size=1) +
    facet_wrap(~Compartment, scales='free_y', ncol=8) +
    ylab(region_to_plot) + 
    theme_bw()->p
    if (region_to_plot == 0){
      p = p + ylab('Illinois')
    }
p
} -> plots

g = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], nrow=12)
ggsave(outfile2, g, width=32, height=24)


### Format for civis ###

covid_set_root(root)

intervention = model_name

agg = plotout %>% select(Date, Region, Compartment, Cases)
agg$intervention=intervention
print('Agg success')
ptilehigh <- function(x){return(quantile(x, 0.975))}
ptilelow <- function(x){return(quantile(x, 0.025))}
ptilemed <- function(x){return(quantile(x, 0.5))}

Rt_all = agg %>% filter(Compartment=='Rt') %>%
    group_by(Date, intervention, Region) %>%
    summarize(rt_median=ptilemed(Cases),
          rt_lower=ptilelow(Cases),
          rt_upper=ptilehigh(Cases))

prev = agg %>% filter(Compartment=='Prevalence') %>%
    group_by(Date, intervention, Region) %>%
    summarize(cases_median=ptilemed(Cases),
          cases_lower=ptilelow(Cases),
          cases_upper=ptilehigh(Cases))
print('Got prevalence')
inc = agg %>% filter(Compartment=='Incidence') %>%
    group_by(Date, intervention, Region) %>%
    summarize(cases_new_median=ptilemed(Cases),
          cases_new_lower=ptilelow(Cases),
          cases_new_upper=ptilehigh(Cases))

deaths = agg %>% filter(Compartment=='nD') %>%
    group_by(Date, intervention, Region) %>%
    summarize(deaths_median=ptilemed(Cases),
           deaths_lower=ptilelow(Cases),
          deaths_upper=ptilehigh(Cases))
    
deaths_det = agg %>% filter(Compartment=='Reported deaths') %>%
    #group_by(SimID, Date, intervention, Region) %>%
    #summarize(Cases = sum(Cases)) %>%
    mutate(Compartment = 'deaths_det') %>%
    #ungroup() %>%
    group_by(Date, intervention, Region) %>%
    summarize(deaths_det_median=ptilemed(Cases),
           deaths_det_lower=ptilelow(Cases),
          deaths_det_upper=ptilehigh(Cases))


hosp = agg %>% filter(Compartment=='IH') %>%
    group_by(Date, intervention, Region) %>%
    summarize(hosp_bed_median=ptilemed(Cases),
          hosp_bed_lower=ptilelow(Cases),
          hosp_bed_upper=ptilehigh(Cases))

icu = agg %>% filter(Compartment=='IC') %>%
    group_by(Date, intervention, Region) %>%
    summarize(icu_median=ptilemed(Cases),
          icu_lower=ptilelow(Cases),
          icu_upper=ptilehigh(Cases))

vent = icu %>% 
    mutate(vent_median=0.74*icu_median,
           vent_lower=0.74*icu_lower,
           vent_upper=0.74*icu_upper) %>%
    select(-icu_median, -icu_lower, -icu_upper)

recov = agg %>% filter(Compartment=='R') %>%
    group_by(Date, intervention, Region) %>%
    summarize(recovered_median=ptilemed(Cases),
          recovered_lower=ptilelow(Cases),
          recovered_upper=ptilehigh(Cases))

final_output = inner_join(prev, inc, by=c('Date', 'intervention', 'Region')) %>% 
    inner_join(deaths, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(deaths_det, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(hosp, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(icu, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(vent,by=c('Date', 'intervention', 'Region')) %>%
    inner_join(recov, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(Rt_all, by=c('Date', 'intervention', 'Region')) %>%
    filter(Date >= as.Date('2020-03-13')) %>% 
    arrange(Region, intervention, Date) %>%
    ungroup()

names(final_output) = c('date', 'scenario_name', 'geography_modeled', names(final_output[4:ncol(final_output)]))
final_output = final_output %>% 
  mutate(geography_modeled=if_else(geography_modeled==0, 'illinois', paste0('covidregion_',geography_modeled))) %>%
  select(date, geography_modeled, scenario_name, names(final_output[4:ncol(final_output)]))
write.csv(final_output, sprintf('%s/uchicago_%s.csv',output_dir, format(Sys.Date(), "%Y%m%d")), row.names=F)
