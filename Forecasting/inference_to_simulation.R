## Source simulation infrastructure
source('simulation_statewide.R')
# CSnippets

rprocess_Csnippet <- read_Csnippet(rprocFile)
rinit_Csnippet <- read_Csnippet(initFile)

# import scaling for nu
nu_scales = read.csv(nu_scales_file)
names(nu_scales) = c('time', 'nu_scale')
row.names(nu_scales) = nu_scales$time

# import fraction underreported covariate
fraction_underreported = read.csv(fraction_underreported_file)
names(fraction_underreported) = c('time', 'frac_underreported')
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
  
  beta_scaling = scales[i]
  
  # Generate scaling for beta
  beta_scales <- get_scale(t_logistic_start = as.numeric(today - reference_date),
                           intervention_lift = as.numeric(intervention_lift_date - reference_date),
                           simstart =  as.numeric(start_projection_date - reference_date),
                           simend = as.numeric(end_projection_date - reference_date),
                           max_scale = beta_scaling
                           )
  
  
  # Define scenario based on beta scale
  if (beta_scaling == 1.0){
    scenario = 100
  } else{
    print(paste0('Scaling is ', beta_scaling))

    scenario = (1 - (beta_scaling - 1)) * 100
    print(paste0('Scenario is ', scenario))
  }
  
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
    pars$nu_m = r[['nu_m']]


    if (pars$n_regions == 4){
        pars$beta2_4 = r[['beta2_4']]
        pars$num_init_4 = r[['num_init_4']]
    }

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
        rinit_Csnippet = rinit_Csnippet
        ) 

    get_obsprob = function(Time){
      pars$nu_3 * pars$theta_test * (1 - fraction_underreported[as.character(Time), 'frac_underreported'] * runif(length(Time), 0.8, 1))
    }

    get_obsprob_nonhosp = function(Time){
        get_fracunder = function(x) {fraction_underreported[as.character(x), 'frac_underreported']}
        get_nu = function(x){nu_scales[as.character(x), 'nu_scale']}
        underreporting = case_when((Time < 75) ~ as.numeric(1 - get_fracunder(Time)),
                                   (Time >= 75) ~ as.numeric(get_nu(Time)))
        pars$theta_test * pars$nu_m * underreporting
    }

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
    print('imposing observation model')
    # Impose observation model on hospitalized deaths
    hd = sims_total %>% filter(Compartment == 'nHD') %>%
      mutate(Cases = rbetabinom(length(Time), Cases, get_obsprob(Time), pars$dispersion))

    # Impose observation model on non-hospitalized deaths
    nh = sims_total %>% filter(Compartment %in% c('nD', 'nHD')) %>%
        spread(Compartment, Cases)  %>%
        mutate(nDnH = nD - nHD) %>% 
        select(-nD, -nHD) %>%
        gather(Compartment, Cases, 'nDnH') %>%
        mutate(Cases = rbetabinom(length(Time), Cases, get_obsprob_nonhosp(Time), pars$dispersion))

    print('Cumulative deaths')
    # Calculate cumulative deaths
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
        filter(Compartment %in% c('nHD', 'nDnH'), Region!="0") %>%
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
  
    sims_total <- rbind(sims_total,
                        hd %>% mutate(Compartment = "Reported hd"),
                        nh %>% mutate(Compartment = "Reported nhd"))

    sims_total$parset = jobid
    sims_total = as.data.frame(sims_total)
    sims_total_combined = rbind(sims_total_combined, sims_total)
  }
  
  
  ## Output files for forecasting hub
  final = final %>% filter(Time <= as.numeric(end_forecast_hub_date - reference_date))
  format_for_covid_hub(final) -> projection_daily
  format_for_covid_hub_week(final) -> projection_weekly
  bind_rows(projection_daily, projection_weekly) -> final_projection
  write.csv(final_projection, paste0(output_path, sprintf("%s-UChicago-CovidIL_%s.csv", today, scenario)))

  ## Output files for plotting
  sims_total_combined$Date = reference_date + sims_total_combined$Time
  write.csv(sims_total_combined %>% filter(Region!='0'), paste0(output_path,sprintf("plotting_output_%s.csv", scenario)))
  write.csv(sims_total_combined %>% filter(Region=='0') %>% select(-Region), paste0(output_path,sprintf("plotting_output_statewide_%s.csv", scenario)))
}
