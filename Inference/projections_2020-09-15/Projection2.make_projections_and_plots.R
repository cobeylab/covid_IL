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

final_date = as.Date('2020-09-01')
root <- '../../'
source(file.path(root, '_covid_root.R'))

covid_set_root(root)
default_par_file = './default_parameter_values.csv'
print(default_par_file)
fit_df = read.csv('./full.final_points.csv')


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
    mutate(confirmed_covid_nonicu = if_else(Date < as.Date('2020-05-06'), NA_integer_, confirmed_covid_nonicu))

R0_values = read.csv(paste0(projection_dir, 'R0_values.csv')) %>%
  mutate(Date=reference_date+time,
    Region=as.character(Region))

## Simulate up to current date

# Loop through and do simulations for each parameter set
for (parset_id in 1:max(fit_df$parset)){
 
    pars_to_use = fit_df %>% 
      filter(parset == parset_id) %>% 
      select(-X) %>% 
      unlist(use.names=T)
    pars[names(pars_to_use)] = pars_to_use
    num_sims = pars$num_sims

    end_projection_date <- as.Date(final_projection_date) #Sys.Date()
    start_projection_date <- simstart + reference_date

    ## Set up simulation timng params
    t0_sim = as.numeric(start_projection_date - reference_date)
    simend = as.numeric(end_projection_date - reference_date)
    pars$tmin=t0_sim
    pars$tmax = simend 

    sim_full = simulate_pomp_covid(
      n_regions = pars$n_regions,
      n_age_groups = n_age_groups,
      nsim=num_sims, 
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
      initialize=T,
      obsnames =observed_names
    )
        sim_full$raw_simulation_output$.id = as.numeric(sim_full$raw_simulation_output$.id)
        alloutputs = process_pomp_covid_output(sim_full, agg_regions=F)
        
        plotout = alloutputs$plotting_output %>% ungroup() %>%    
            mutate(Date=reference_date + Time)
        df_infections_summary <- plotout %>%
            group_by(Time, Compartment, Region) %>%
            filter(Compartment %in% c('E', 'A', 'P', 'IS', 'IM', 'IH', 'IC')) %>%
            group_by(SimID, Date, Region) %>%
            arrange(Time) %>%
            summarise(Compartment = 'Prevalence',
                      Cases = sum(Cases)) %>%
          ungroup()
        
        plotout = bind_rows(plotout, df_infections_summary) %>% 
            select(-Time)
        
        statewide = plotout %>%
            group_by(Date, SimID, Compartment) %>%
            summarize(Cases=sum(Cases)) %>%
            mutate(Region='0') %>%
            ungroup()
        
        plotout = bind_rows(plotout, statewide) %>% filter(!is.na(Compartment))
        plotout$parset = parset_id

        Rtdf = add_R0_to_output(plotout, R0_values, parset)
        plotout = rbind(plotout, Rtdf)
        
        
        write.csv(plotout, paste0(model_name, '_parset_', parset_id, '.csv'), row.names=F)
        rm(sim_full)
        
    }

file_list <- list.files(pattern = paste0(model_name, '_parset_.*.csv'))
cons <- rbindlist(lapply(file_list, fread))
plotout <- as.data.frame(cons) %>% mutate(Date=as.Date(Date))
plotout = plotout %>% mutate(covid_region = Region)
write.csv(plotout, paste0(output_dir,model_name, '_all.csv'))