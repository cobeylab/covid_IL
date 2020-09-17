library(tidyverse)
library(lubridate)



civis_template <- as_tibble(read.csv('capacity_weekday_average_20200915.csv')) %>%
  mutate(date_window_upper_bound = ymd(date_window_upper_bound))
start_date <- ymd(unique(civis_template$date_capacity_run)) # Does hospital capacitity exceed threshold at some point beginning on this date?
stopifnot(length(start_date) == 1)

print('Reading projections...')
simulations <- as_tibble(read.csv(sprintf('baseline_all.consolidated_%s.csv', format(Sys.Date(), "%Y%m%d")))) %>%
  mutate(Date = ymd(Date),
        set = 1)

print('Calculating % simulations exceeding capacity')
get_fraction_sims_exceeding_capacity <- function(date_window_upper_bound, start_date, simulations, civis_template){
  
  thresholds <- civis_template %>% filter(date_window_upper_bound == !!date_window_upper_bound) %>%
    select(date_window_upper_bound, geography_modeled, resource_type, avg_resource_available_prev2weeks, overflow_threshold_percent) %>%
    separate(geography_modeled, into = c('useless_prefix','covid_region'), sep = '_') %>%
    select(-useless_prefix) %>%
    mutate(covid_region = as.integer(covid_region))
  
  occupancy <- simulations %>%
    filter(covid_region != 0) %>%
    filter(Date >= start_date, Date <= date_window_upper_bound, Compartment %in% c('IH','IC')) %>%
    group_by(SimID, parset, set, Region, Date, Compartment) %>%
    summarise(occupancy = sum(Cases)) %>%
    mutate(resource_type = case_when(
      Compartment == 'IC' ~ 'icu_availforcovid',
      Compartment == 'IH' ~ 'non_icu_only'
    )) %>% ungroup()
  
  # Add total beds as another resource type, keep just total and ICU for merging with civis template
  total_occupancy <- occupancy %>% group_by(SimID, parset, set, Region, Date) %>%
    summarise(occupancy = sum(occupancy)) %>% mutate(resource_type = 'hb_availforcovid') %>%
    ungroup()
  
  
  occupancy <- bind_rows(occupancy, total_occupancy) %>%
    arrange(SimID, parset, set, Region, Date) %>%
    filter(resource_type %in% c('icu_availforcovid','hb_availforcovid')) %>% 
    select(-Compartment)
  
  # For each simulation, find maximum occupancy between start_date and date_window_upper_bound
  max_occupancy <- occupancy %>% group_by(SimID, parset, set, Region, resource_type) %>%
    summarise(max_occupancy = max(occupancy)) %>%
    rename(covid_region = Region) %>% 
    ungroup()
  
  # For each region, count how many simulations exceed the 75% and 100% thresholds
  threshold_check <- full_join(max_occupancy, 
            thresholds %>% select(covid_region, resource_type, avg_resource_available_prev2weeks, overflow_threshold_percent),
            by = c('covid_region', 'resource_type')) %>%
    group_by(covid_region, resource_type, overflow_threshold_percent) %>%
    summarise(percent_of_simulations_that_exceed = sum(max_occupancy > avg_resource_available_prev2weeks)/length(max_occupancy)) %>%
    ungroup()
  
  # As a test, check that the number of simulations exceeding the 100% threshold is smaller than or equal to the number exceeding 75%
  makes_sense <- threshold_check %>% group_by(covid_region, resource_type) %>%
    summarise(makes_sense = percent_of_simulations_that_exceed[overflow_threshold_percent == 1] <= percent_of_simulations_that_exceed[overflow_threshold_percent == 0.75]) %>%
    pull(makes_sense)
  stopifnot(makes_sense)
  
  return(left_join(thresholds, threshold_check, by = c('covid_region','resource_type','overflow_threshold_percent')))
  
}

threshold_checks <-  lapply(as.list(unique(civis_template$date_window_upper_bound)), get_fraction_sims_exceeding_capacity,
       start_date = start_date, simulations = simulations, civis_template = civis_template)

threshold_checks <- bind_rows(threshold_checks) %>%
  mutate(covid_region = paste0('covidregion_', covid_region)) %>%
  rename(geography_modeled = covid_region)

filled_out_template <- left_join(civis_template %>% select(-percent_of_simulations_that_exceed),
          threshold_checks,
          by = c('date_window_upper_bound','geography_modeled','resource_type','avg_resource_available_prev2weeks',
                 'overflow_threshold_percent')) %>%
  mutate(scenario_name = 'baseline')

write.csv(filled_out_template, sprintf('uchicago_hospitaloverflow_%s.csv',  format(Sys.Date(), "%Y%m%d")), row.names = F)
print('Done.')
