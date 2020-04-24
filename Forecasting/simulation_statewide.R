#' Read a whole file.
read_text_file <- function(path) {
  readChar(path, file.info(path)$size)
}


#' Read a Csnippet
read_Csnippet <- function(path) {
  library(pomp)
  Csnippet(read_text_file(path))
}

#' Simulate COVID pomp model
simulate_pomp_covid <- function(
  n_age_groups,
  n_regions,
  nsim, 
  input_params, 
  delta_t, 
  contacts, 
  population_list, 
  nu_scales, 
  frac_underreported,
  seed = NULL,
  format = 'data.frame',
  obsnames = NULL,
  rprocess_Csnippet= NULL,
  initialize = T,
  rinit_Csnippet= NULL,
  pfilter_traj_params = NULL
) {
  library(pomp)
  
  # Initialize states
  subcompartment_df <- simulate_pomp_covid__init_subcompartment_df(input_params)

  # Intervention df
  intervention_df <- simulate_pomp_covid__init_intervention_df(input_params)
  n_interventions <- nrow(intervention_df)

  # Set up parameters for pomp
  params <- simulate_pomp_covid__init_parameters(
    n_age_groups, n_interventions, input_params, contacts, population_list
  )

  # Set up covariate table for pomp
  covar_table_interventions <- simulate_pomp_covid__init_covariate_table(
    input_params, intervention_df, nu_scales, frac_underreported
  )
  covar_table <- covariate_table(
    covar_table_interventions,
    order = "constant",
    times = "time"
  )
  # Actually call pomp
  state_names <- simulate_pomp_covid__init_state_names(n_age_groups, n_regions, subcompartment_df)
  
  ## Set up accumulator variables by age group and then by regions
  acc_mild_age <- sprintf("new_mild_infections_%d",c(1:n_age_groups))
  acc_new_severe <- sprintf("new_severe_infections_%d",c(1:n_age_groups))
  acc_nD <- sprintf("new_deaths_%d",c(1:n_age_groups))
  acc_IH1 <- sprintf("new_IH1_%d",c(1:n_age_groups))
  acc_IC2 <- sprintf("new_IC2_%d",c(1:n_age_groups))
  acc_IC3 <- sprintf("new_IC3_%d",c(1:n_age_groups))
  acc_IH4 <- sprintf("new_IH4_%d",c(1:n_age_groups))
  acc_INC <- sprintf("Inc_%d", c(1:n_age_groups))
  
  accum_names <- array()
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_mild_age, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_new_severe, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_nD, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_IH1, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_IC2, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_IC3, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_IH4, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_INC, "_", i))
  }
  
  accum_names <- accum_names[!is.na(accum_names)]
  times <- c(ceiling(input_params$tmin):ceiling(input_params$tmax)) 
  
  if(initialize == T){
  output_df <- simulate(
    nsim = nsim,
    seed = seed,
    times = times,
    t0 = times[1],
    rinit = rinit_Csnippet,
    rprocess = euler(rprocess_Csnippet, delta.t = delta_t),
    #rmeasure = rmeasure_Csnippet,
    
    params = params,
    covar = covar_table,
    
    statenames = c(state_names, accum_names),
    paramnames = names(params),
    accumvars = accum_names,
    #obsnames = c(
     # paste0("ObsCases_", c(1:n_age_groups)), paste0("ObsDeaths_",c(1:n_age_groups))
    #),
    
    format = format
  )
  }
  if(initialize == F){
    output_df <- simulate(
      nsim = 1,
      seed = seed,
      times = times,
      t0 = times[1],
      #rinit = rinit_Csnippet,
      rprocess = euler(rprocess_Csnippet, delta.t = delta_t),
      rmeasure = rmeasure_Csnippet,
      
      params = pfilter_traj_params,
      covar = covar_table,
      
      statenames = c(state_names, accum_names),
      paramnames = names(params),
      accumvars = accum_names,
      obsnames = c(
        paste0("ObsCases_", c(1:n_age_groups)), paste0("ObsDeaths_",c(1:n_age_groups))
      ),
      
      format = format
    )
  }
  list(
    n_age_groups = n_age_groups,
    raw_simulation_output = output_df,
    params = params,
    state_names = state_names,
    interventions = covar_table_interventions
  )
}

simulate_pomp_covid__init_subcompartment_df <- function(input_params) {
  with(input_params, {
    data.frame(
      state = c("E","A","P","IM","IS", "IH1_", "IH2_", "IH3_", "IC2_", "IC3_", "IH4_"),
      index = c(
        alpha_E,
        alpha_A, alpha_P,
        alpha_IM, alpha_IS,
        alpha_IH1, alpha_IH2, alpha_IH3,
        alpha_IC2, alpha_IC3,
        alpha_IH4
      )
    )
  })
}

simulate_pomp_covid__init_state_names <- function(n_age_groups, n_regions, subcompartment_df) {
  ## State variables with variable classes per age group
  state_names_S_age = c(sprintf("S_%d",c(1:n_age_groups)))
  state_names_R_age = c(sprintf("R_%d",c(1:n_age_groups)))
  state_names_D_age <- c(sprintf("D_%d",c(1:n_age_groups)))
  
  state_names_S <- array()
  state_names_R <- array()
  state_names_D <- array()
  for(i in c(1:n_regions)){
    state_names_S <- c(state_names_S, paste0(state_names_S_age, "_", i))
    state_names_R <- c(state_names_R, paste0(state_names_R_age, "_", i))
    state_names_D <- c(state_names_D, paste0(state_names_D_age, "_", i))
  }
    
  state_names <- c(state_names_S, state_names_R, state_names_D)

  for(i in c(1:nrow(subcompartment_df))) {
      state_var = as.character(subcompartment_df[i,]$state)
      n_subcompartments = subcompartment_df[i,]$index
      state_names_var_age <- array()
      for(k in c(1:n_age_groups)){
        for(j in c(1:n_subcompartments)) {
          state_names_var_age = c(state_names_var_age, paste0(state_var,j,"_",k))
        }
      }
      state_names_var_age <- state_names_var_age[!is.na(state_names_var_age)]
      state_names_var_age_region <- array()
      for(z in c(1:n_regions)){
         state_names_var_age_region <- c(state_names_var_age_region, paste0(state_names_var_age,"_",z))
      }
      state_names_var_age_region <- state_names_var_age_region[!is.na(state_names_var_age_region)]
      state_names <- c(state_names, state_names_var_age_region)
  }

  state_names <- state_names[!is.na(state_names)]
  state_names
}


simulate_pomp_covid__init_parameters <- function(
  n_age_groups, n_interventions, input_params, contacts, population_list
) {
  with(input_params, {
    # Set up parameters for model
    params = c(
      "n_age_groups" = n_age_groups,
      "n_regions" = n_regions,
      "n_interventions" = n_interventions,
      "inference" = inference,
      "simulation" = simulation,
      "sigma" = 1/inv_sigma,
      "alpha_E"= alpha_E,
      "alpha_P"= alpha_P,
      "alpha_A"= alpha_A,
      "alpha_IM"= alpha_IM,
      "alpha_IS"= alpha_IS,
      "alpha_IH1" = alpha_IH1, 
      "alpha_IH2" = alpha_IH2,
      "alpha_IH3" = alpha_IH3, 
      "alpha_IC2" = alpha_IC2,
      "alpha_IC3" = alpha_IC3,
      "alpha_IH4" = alpha_IH4,
      "beta1" = beta1,
      #"beta2" = beta2,
      "theta_test"=theta_test,
      #"num_init"=num_init,
      "nu_1"= nu_1,
      "nu_3"= nu_3,
      "nu_m"= nu_m,
      "dispersion" = dispersion,
      "frac_severe_init" = frac_severe_init,
      "t_reporting_adjustment" = t_reporting_adjustment,
      "lower_bound_reporting_uncertainty" = lower_bound_reporting_uncertainty
    )
    
    for(i in c(1:length(population_list))){
      params[paste0(paste0("N_",c(1:n_age_groups)),"_",i)] = population_list[[i]]$POPULATION
    }
    
    for(i in c(1:length(population_list))){
      params[paste0(paste0("age_dist_",c(1:n_age_groups)),"_",i)] = population_list[[i]]$POPULATION / sum(population_list[[i]]$POPULATION)
    }
    
    # Beta 2's by region
    params[paste0("beta2_",c(1:n_regions))] = c(beta2_1, beta2_2, beta2_3)
    
    # Num init by region 
    params[paste0("num_init_",c(1:n_regions))] = c(num_init_1, num_init_2, num_init_3)
    
    # Rates
    params[paste0("eta_",c(1:n_age_groups))] = rep(1/inv_eta, n_age_groups)
    params[paste0("zeta_s_",c(1:n_age_groups))] = rep(1/inv_zeta_s, n_age_groups)
    params[paste0("zeta_h_",c(1:n_age_groups))] = rep(1/inv_zeta_h, n_age_groups)
    params[paste0("zeta_c_",c(1:n_age_groups))] = rep(1/inv_zeta_c, n_age_groups)
    params[paste0("gamma_m_",c(1:n_age_groups))] = rep(1/inv_gamma_m, n_age_groups)
    params[paste0("gamma_c_",c(1:n_age_groups))] = rep(1/inv_gamma_c, n_age_groups)
    params[paste0("gamma_h_",c(1:n_age_groups))] = rep(1/inv_gamma_h, n_age_groups)
    params[paste0("mu_c_",c(1:n_age_groups))] = rep(1/inv_mu_c, n_age_groups)
    params[paste0("mu_h_",c(1:n_age_groups))] = rep(1/inv_mu_h, n_age_groups)
    
    # Probabilities
    params[paste0("rho_",c(1:n_age_groups))] = unlist(
      input_params[paste0('rho_', c(1:n_age_groups))]
    )
    params[paste0("kappa_",c(1:n_age_groups))] = unlist(
      input_params[paste0('kappa_', c(1:n_age_groups))]
    )
    params[paste0("psi1_",c(1:n_age_groups))] = unlist(
      input_params[paste0('psi1_', c(1:n_age_groups))]
    )
    params[paste0("psi2_",c(1:n_age_groups))] = unlist(
      input_params[paste0('psi2_', c(1:n_age_groups))]
    )
    params[paste0("psi3_",c(1:n_age_groups))] = unlist(
      input_params[paste0('psi3_', c(1:n_age_groups))]
    )
    #params[paste0("nu_s_",c(1:n_age_groups))] = rep(nu_s, n_age_groups)
    #params[paste0("nu_m_",c(1:n_age_groups))] = rep(nu_m, n_age_groups)
  
    for(i in c(1:n_age_groups)){
      for(j in c(1:n_age_groups)){
        params[paste0("C_home_",j,"_",i)] = 0
      }
    }
    for(i in c(1:n_age_groups)){
      for(j in c(1:n_age_groups)){
        params[paste0("C_work_",j,"_",i)] = 0
      }
    }
    
    for(i in c(1:n_age_groups)){
      for(j in c(1:n_age_groups)){
        params[paste0("C_school_",j,"_",i)] = 0
      }
    }
    
    for(i in c(1:n_age_groups)){
      for(j in c(1:n_age_groups)){
        params[paste0("C_other_",j,"_",i)] = 0
      }
    }
    params[grep("C_home_",names(params))] = contacts$home
    params[grep("C_work_",names(params))] = contacts$work
    params[grep("C_school_",names(params))] = contacts$school
    params[grep("C_other_",names(params))] = contacts$other
    
    params
  })
}

simulate_pomp_covid__init_intervention_df <- function(input_params) {
  # Generate intervention data frame 
  with(input_params, {
    intervention_df <- data.frame(
      t_start = c(t_start), # example with 2 interventions
      t_end = c(t_end),
      scale_work = c(scalework), # scale this type of contacts
      scale_home = c(scalehome),
      scale_school = c(scaleschool),
      scale_other = c(scaleother)
    )
  })
}


simulate_pomp_covid__init_covariate_table <- function(input_params, intervention_df, nu_scales, frac_underreported) {
  n_interventions <- nrow(intervention_df)
  
  # Specify interventions via covariate table 
  name_vec <- c(
    paste0("t_start_intervention_",c(1:n_interventions)), 
    paste0("t_end_intervention_",c(1:n_interventions)),
    paste0("scale_work_intervention_",c(1:n_interventions)),
    paste0("scale_home_intervention_",c(1:n_interventions)),
    paste0("scale_school_intervention_",c(1:n_interventions)),
    paste0("scale_other_intervention_",c(1:n_interventions))
  )
  
  intervention_matrix <- matrix(
    rep(c(
      intervention_df$t_start, 
      intervention_df$t_end, 
      intervention_df$scale_work,
      intervention_df$scale_home,
      intervention_df$scale_school,
      intervention_df$scale_other
    ), input_params$tmax),
    ncol = length(name_vec),
    byrow = T
  )
  
  covar_table_interventions <- data.frame(cbind(
    1:input_params$tmax,
    matrix(
      rep(c(
        intervention_df$t_start, 
        intervention_df$t_end, 
        intervention_df$scale_work,
        intervention_df$scale_home,
        intervention_df$scale_school,
        intervention_df$scale_other
      ), input_params$tmax),
      ncol = length(name_vec),
      byrow = T)
  ))
  names(covar_table_interventions) <- c("time", name_vec)
  
  covar_table_interventions$nu_scale = nu_scales[covar_table_interventions$time, 'nu_scale']
  
  covar_table_interventions$frac_underreported = frac_underreported[covar_table_interventions$time, 'frac_underreported']

  covar_table_interventions
}

#' Make output table with the following columns:
#' Simulation ID, Time, Compartment, Age group, Cases
 
process_pomp_covid_output <- function(sim_result) {
  select <- dplyr::select
  rename <- dplyr::rename
  summarize <- dplyr::summarise
  contains <- dplyr::contains
  
  n_age_groups <- sim_result$n_age_groups
  df_sim <- sim_result$raw_simulation_output
  params <- sim_result$params

  df_sim_output <- df_sim %>% 
      melt(id.vars = c(".id", "time")) %>% 
      rename(
          Compartment = variable,
          Cases = value
      )
  df_sim_output %>% mutate(Compartment = case_when(
      startsWith(as.character(Compartment), "S") ~ "S",
      startsWith(as.character(Compartment), "E") ~ "E",
      startsWith(as.character(Compartment), "P") ~ "P",
      startsWith(as.character(Compartment), "A") ~ "A",
      startsWith(as.character(Compartment), "IS") ~ "IS",
      startsWith(as.character(Compartment), "IM") ~ "IM",
      startsWith(as.character(Compartment), "IH1") ~ "IH",
      startsWith(as.character(Compartment), "IH2") ~ "IH",
      startsWith(as.character(Compartment), "IH3") ~ "IH",
      startsWith(as.character(Compartment), "IH4") ~ "IH",
      startsWith(as.character(Compartment), "IC2") ~ "IC",
      startsWith(as.character(Compartment), "IC3") ~ "IC",
      startsWith(as.character(Compartment), "R") ~ "R",
      startsWith(as.character(Compartment), "D") ~ "D",
      startsWith(as.character(Compartment), "ObsCases") ~ "O",
      startsWith(as.character(Compartment), "ObsDeaths") ~ "OD",
      startsWith(as.character(Compartment), "Inc") ~ "Incidence",
      startsWith(as.character(Compartment), "new_deaths") ~ "nD",
      startsWith(as.character(Compartment), "new_mild") ~ "nM",
      startsWith(as.character(Compartment), "new_severe") ~ "nS"
  )) -> df_sim_output
  names(df_sim_output) = c('SimID', 'Time', 'Compartment', 'Cases')
  df_sim_output %>%
    group_by(SimID, Time, Compartment) %>%
    summarize(Cases = sum(Cases)) -> df_sim_output
  c(
    sim_result,
    list(plotting_output = df_sim_output)
  )
}

