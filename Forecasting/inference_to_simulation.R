## Source simulation infrastructure
source('simulation_statewide.R')
# CSnippets

rprocess_Csnippet <- read_Csnippet(rprocFile)
rinit_Csnippet <- read_Csnippet(initFile)
rmeasure_Csnippet <- read_Csnippet(rmeasFile)

# import scaling for nu
nu_scales = read.csv(nu_scales_file)
names(nu_scales) = c('time', 'nu_scale')
row.names(nu_scales) = nu_scales$time

# import fraction underreported covariate
fraction_underreported = read.csv(fraction_underreported_file)
row.names(fraction_underreported) = fraction_underreported$time

# Import contact matrix data
load(contact_filename)

print('Reading in population files')
# Make population a vector that has the regions concatenated together
population1 = read.csv(population_filename_1)
colnames(population1) <- c("AGE_GROUP_MIN", "POPULATION")

population2 = read.csv(population_filename_2)
colnames(population2) <- c("AGE_GROUP_MIN", "POPULATION")

population3 = read.csv(population_filename_3)
colnames(population3) <- c("AGE_GROUP_MIN", "POPULATION")

population4 = read.csv(population_filename_4)
colnames(population4) <- c("AGE_GROUP_MIN", "POPULATION")

n_age_groups = nrow(population1)

### User-specification of parameters and interventions for model 
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)

## Set up simulation timng params
t0_sim = as.numeric(start_projection_date - reference_date)
simend = as.numeric(end_projection_date - reference_date)
pars$tmin=t0_sim
pars$tmax = simend 

pars <- add_interventions(intervention_table_filename, pars)


# ASSUME FOLLOWING ORDER: NORTH-CENTRAL, CENTRAL, NORTHEAST, SOUTHERN 
if (pars$n_regions == 3){
    population2$POPULATION = population2$POPULATION + population4$POPULATION

    population_list=list(population1, population2, population3)
} else if (pars$n_regions == 4){
    population_list=list(population1, population2, population3, population4)
}


for(i in c(1:length(scales))){
  
  beta_scaling = as.numeric(scales[[i]])
  
  # Generate scaling for beta
  beta_scales <- get_scale(t_logistic_start = as.numeric(today - reference_date),
                           intervention_lift = as.numeric(intervention_lift_date - reference_date),
                           simstart =  as.numeric(start_projection_date - reference_date),
                           simend = as.numeric(end_projection_date - reference_date),
                           max_scales = beta_scaling
                           )
  
  
  # Define scenario based on beta scale
  scenario = scenarios[i]
  
  cat("running scenario ", scenario, "\n")
  temp_scales = beta_scales
  
  pars$beta_noise_amplitude = beta_noise_amplitude
  
  design <- read.csv(input_points_file) %>% 
    arrange(-loglik) 
  
  final = data.frame()
  sims_total_combined = data.frame()
  
  for(jobid in c(1:nrow(design))){
    r = design[jobid, ]
    pars$beta1 = r[['beta1']]
    pars$beta2_1 = r[['beta2_1']]
    pars$beta2_2 = r[['beta2_2']]
    pars$beta2_3 = r[['beta2_3']]
    pars$num_init_1 = r[['num_init_1']]
    pars$num_init_2 = r[['num_init_2']]
    pars$num_init_3 = r[['num_init_3']]
    pars$region_non_hosp_1 = r[['region_non_hosp_1']]
    pars$region_non_hosp_2 = r[['region_non_hosp_2']]
    pars$region_non_hosp_3 = r[['region_non_hosp_3']]
    pars$region_non_hosp_4 = r[['region_non_hosp_4']]
    pars$beta2_4 = r[['beta2_4']]
    pars$num_init_4 = r[['num_init_4']]


    n_projections = r[['num_sims']]
    print(deltaT)
    sim_full <- simulate_pomp_covid(
        n_age_groups=nrow(population_list[[1]]),
        n_regions = pars$n_regions,
        nsim= n_projections, 
        input_params=pars,
        delta_t=deltaT,
        contacts=pomp_contacts,
        population_list = population_list,
        nu_scales=nu_scales,
        beta_scales= beta_scales,
        frac_underreported = fraction_underreported,
        rprocess_Csnippet = rprocess_Csnippet,
        initialize = T,
        rinit_Csnippet = rinit_Csnippet,
        rmeasure_Csnippet=rmeasure_Csnippet,
        obsnames=c(paste0('ObsHospDeaths_', 1:4), paste0('ObsNonHospDeaths_', 1:4), paste0('ObsICU_', 1:4)) 
        ) 

    ## Generate the direct simulation output
    sim_full %>% process_pomp_covid_output(agg_regions = regional_aggregation)  -> simout

    sims_total = as.data.frame(simout$plotting_output) %>% 
      filter(!is.na(Compartment))

    statewide = sims_total %>% 
        group_by(SimID, Time, Compartment) %>% 
        summarize(Cases=sum(Cases)) %>% 
        ungroup() %>% 
        mutate(Region='0')

    sims_total = bind_rows(sims_total, statewide)

    print(head(sims_total))

    #print('Cumulative deaths')
    # Calculate cumulative deaths
    nh = sims_total %>% filter(Compartment=='Reported nhd')
    hd = sims_total %>% filter(Compartment=='Reported hd')
    
    cumulative_nh = nh %>% 
      group_by(SimID, Region) %>%
      arrange(Time) %>%
      mutate(Compartment='DnH',
             Cases=cumsum(Cases)) %>%
      ungroup()
    cumulative_hd = hd %>% 
      group_by(SimID, Region) %>%
      arrange(Time) %>%
      mutate(Compartment='D',
             Cases=cumsum(Cases)) %>%
      ungroup()

    # For forecast hub
      new_deaths <- bind_rows(hd, nh) %>%
        filter(Compartment %in% c('Reported hd', 'Reported nhd'), Region!="0") %>%
        group_by(SimID, Time) %>%
        summarize(Cases=sum(Cases)) %>%
        mutate(Compartment='nD') %>%
        ungroup()
      tot_deaths <- bind_rows(cumulative_nh, cumulative_hd) %>%
        filter(Compartment %in% c('D', 'DnH'), Region!="0") %>%
        group_by(SimID, Time) %>%
        summarize(Cases=sum(Cases)) %>%
        mutate(Compartment='D') %>%
        ungroup()
    sims <- bind_rows(new_deaths, tot_deaths)
    sims$parset = jobid
    sims = as.data.frame(sims)
    final <- rbind(final, sims)
  

    sims_total$parset = jobid
    sims_total = as.data.frame(sims_total)
    sims_total_combined = rbind(sims_total_combined, sims_total)
  }
  
  
  ## Output files for forecasting hub
  #final = final %>% filter(Time <= as.numeric(end_forecast_hub_date - reference_date))
  #format_for_covid_hub(final) -> projection_daily
  #format_for_covid_hub_week(final) -> projection_weekly
  #bind_rows(projection_daily, projection_weekly) -> final_projection
  #write.csv(final_projection, paste0(output_path, sprintf("%s-UChicago-CovidIL_%s.csv", today, scenario)))

  ## Output files for plotting
  sims_total_combined$Date = reference_date + sims_total_combined$Time
  write.csv(sims_total_combined %>% filter(Region!='0'), paste0(output_path,sprintf("plotting_output_%s.csv", scenario)))
  write.csv(sims_total_combined %>% filter(Region=='0') %>% select(-Region), paste0(output_path,sprintf("plotting_output_statewide_%s.csv", scenario)))
}
