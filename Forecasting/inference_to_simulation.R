## Source simulation infrastructure
source('simulation_statewide.R')

intervention_start = as.numeric(intervention_start_date - reference_date)
scale_beta_start = as.numeric(scale_beta_start_date - reference_date)

# CSnippets

rprocess_Csnippet <- read_Csnippet(rprocFile)
rinit_Csnippet <- read_Csnippet(initFile)

# import scaling for nu
nu_scales = read.csv(nu_scales_file)
row.names(nu_scales) = nu_scales$time

# generate scaling for beta
beta_scales <- data.frame(time = c(1:as.numeric(end_projection_date - reference_date)),
                          scale_beta = NA
) 
beta_scales[beta_scales$time < scale_beta_start,]$scale_beta = 1.0
beta_scales[beta_scales$time >= scale_beta_start,]$scale_beta = beta_scaling

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

### User-specification of parameters and interventions for model 
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)

## Set intervention start
pars$t_start = intervention_start

## Set intervention strength

  pars$scalework =scalework
  pars$scaleschool=scaleschool
  pars$scaleother=scaleother
  pars$scalehome=scalehome


t0_sim = as.numeric(start_projection_date - reference_date)
pars$tmin=t0_sim
pars$tmax <- as.numeric(end_projection_date - reference_date)
pars$t_end <- as.numeric(intervention_end_date - reference_date)
input_points <- read.csv(input_points_file) %>% 
  arrange(-loglik) 
  

for (row in c(1:nrow(input_points))){
  cat("row is ", row, "\n")
  r = input_points[row, ]
  pars$beta1 = r[['beta1']]
  pars$beta2_1 = r[['beta2_1']]
  pars$beta2_2 = r[['beta2_2']]
  pars$beta2_3 = r[['beta2_3']]
  pars$num_init_1 = r[['num_init_1']]
  pars$num_init_2 = r[['num_init_2']]
  pars$num_init_3 = r[['num_init_3']]
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
      beta_scales= beta_scales,
      frac_underreported = fraction_underreported,
      rprocess_Csnippet = rprocess_Csnippet,
      initialize = T,
      rinit_Csnippet = rinit_Csnippet
      ) 
  
  ## Generate the direct simulation output
  sim_full %>% process_pomp_covid_output(agg_regions = regional_aggregation)  -> d
  
  d_plot = as.data.frame(d$plotting_output) %>% filter(!is.na(Compartment))
  
  ## Generate new, lagged deaths that occurred outside of the hospital ---------------------------------------------------------------------
  df_new_nonhosp_deaths <- add_non_hospitalized_deaths(sim_full = sim_full,
                                                       pars = pars,
                                                       regional_aggregation = regional_aggregation) %>% 
    as.data.frame()
  
  df_new_nonhosp_deaths_leading_zeros <- rbind(df_new_nonhosp_deaths,
                                 data.frame(expand.grid(Time = as.numeric(c(min(sim_full$raw_simulation_output$time):(min(sim_full$raw_simulation_output$time) + time_lag_non_hosp_deaths-1))),
                                                        SimID = unique(df_new_nonhosp_deaths$SimID),
                                                        Compartment = as.character(unique(df_new_nonhosp_deaths$Compartment)),
                                                        Cases = 0))) %>% arrange(Time)

  d_plot <- rbind(d_plot, df_new_nonhosp_deaths_leading_zeros)

  d_plot$Date = as.Date("2020-01-14") + d_plot$Time
  d_plot$Parset = row
    
  filename = paste0(output_filename,"_", row, ".csv")
  write.csv(d_plot, file=filename)
  
}


