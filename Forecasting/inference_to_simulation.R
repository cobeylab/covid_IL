## Source simulation infrastructure
source(covid_get_path('POMP/simulation/statewide/simulation_statewide.R'))
source(covid_get_path('POMP/inference/statewide/inference_statewide.R'))

# CSnippets

rprocess_Csnippet <- read_Csnippet(rprocFile)
rinit_Csnippet <- read_Csnippet(initFile)

# import scaling for nu
nu_scales = read.csv(nu_scales_file)
row.names(nu_scales) = nu_scales$time

# import fraction underreported covariate
fraction_underreported = read.csv(fraction_underreported_file)
row.names(fraction_underreported) = fraction_underreported$time

# Import contact matrix data
load(contact_filename)

# Make population a vector that has the regions concatenated together
population1 = read.csv(population_filename_1)
colnames(population1) <- c("AGE_GROUP_MIN", "POPULATION")

population2 = read.csv(population_filename_2)
colnames(population2) <- c("AGE_GROUP_MIN", "POPULATION")

population3 = read.csv(population_filename_3)
colnames(population3) <- c("AGE_GROUP_MIN", "POPULATION")

population_list=list(population1, population2, population3)

n_age_groups = nrow(population1)

intervention_start = as.numeric(intervention_start_date - reference_date)


### User-specification of parameters and interventions for model 
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)

## Set intervention start
pars$t_start = intervention_start

## Set intervention strength
if(modify_interventions){
  pars$scalework =scalework
  pars$scaleschool=scaleschool
  pars$scaleother=scaleother
  pars$scalehome=scalehome
}

pars$tmin=t0_sim
pars$tmax <- as.numeric(end_projection_date - reference_date)
pars$t_end <- as.numeric(intervention_end_date - reference_date)
input_points <- read.csv(input_points_file)

for (row in c(1:nrow(input_points))){
  cat("row is ", row, "\n")
  r = input_points[row, ]
  pars$beta1 = r[['beta1']]
  pars$beta2_1 = r[['beta2_1']]
  pars$beta2_2 = r[['beta2_2']]
  pars$beta2_3 = r[['beta2_2']]
  pars$num_init_1 = r[['num_init_1']]
  pars$num_init_2 = r[['num_init_2']]
  pars$num_init_1 = r[['num_init_3']]
  n_projections = r[['num_sims']]
  
  sim_full <- simulate_pomp_covid(
      n_age_groups=nrow(population_list[[1]]),
      n_regions = 3,
      nsim= n_projections, 
      input_params=pars,
      delta_t=deltaT,
      contacts=pomp_contacts,
      population_list = population_list,
      nu_scales=nu_scales,
      frac_underreported = fraction_underreported,
      rprocess_Csnippet = rprocess_Csnippet,
      initialize = T,
      rinit_Csnippet = rinit_Csnippet
      ) 
  
  ## Generate the direct simulation output
  sim_full %>% process_pomp_covid_output -> d
  
  d_plot = as.data.frame(d$plotting_output) %>% filter(!is.na(Compartment))
  
  
  ## Generate new, lagged deaths that occurred outside of the hospital ---------------------------------------------------------------------
  sim_raw_IM <- sim_full$raw_simulation_output %>% select(contains(paste0("IM", pars$alpha_IM)))
  
  kappa_vec <- pars[grep("kappa",names(pars))] %>% unlist()
  psi1_vec <- pars[grep("psi1",names(pars))] %>% unlist()
  psi2_vec <- pars[grep("psi2",names(pars))] %>% unlist()
  psi3_vec <- (1-psi1_vec - psi2_vec)
  names(psi3_vec) <- paste0("psi3_", c(1:n_age_groups))
  psi_non_hosp_vec <- (kappa_vec*psi3_vec*frac_non_hosp)/(1-frac_non_hosp)*(1-kappa_vec)
  names(psi_non_hosp_vec) <- paste0("psi_non_hosp_", c(1:n_age_groups))
  
  df_new_nonhosp_deaths <- data.frame(Time =sim_full$raw_simulation_output$time + time_lag_non_hosp_deaths) %>% 
    mutate(SimID = sim_full$raw_simulation_output$.id)
  for(i in c(1:n_age_groups)){
    df_new_nonhosp_deaths[,paste0("NHD_",i)] <- round(sim_raw_IM[,i]*psi_non_hosp_vec[i])
  }
  
  df_new_nonhosp_deaths <- df_new_nonhosp_deaths %>% 
    melt(id.vars = c("SimID", "Time")) %>% 
    rename(
      Compartment = variable,
      Cases = value
    )
  
    df_new_nonhosp_deaths <-  df_new_nonhosp_deaths %>% mutate(Compartment = case_when(
    startsWith(as.character(Compartment), "NHD") ~ "nDnH")) %>%  
    group_by(SimID, Time, Compartment) %>%
    filter(!is.na(Compartment)) %>% 
    summarize(Cases = sum(Cases)) %>% 
      as.data.frame()
    
    
    df_new_nonhosp_deaths <- rbind(df_new_nonhosp_deaths,
                                   expand.grid(Time = c(min(sim_full$raw_simulation_output$time):(min(sim_full$raw_simulation_output$time) + time_lag_non_hosp_deaths-1)),
                                               SimID = unique(df_new_nonhosp_deaths$SimID),
                                               Compartment = unique(df_new_nonhosp_deaths$Compartment),
                                               Cases = 0)) %>% 
      arrange(Time,SimID)

    
    d_plot <- rbind(d_plot, df_new_nonhosp_deaths)
    
    d_plot$Date = as.Date("2020-01-14") + d_plot$Time
    d_plot$Parset = row
    
    filename = paste0(output_filename,"_", row, ".csv")
    write.csv(d_plot, file=filename)
}


