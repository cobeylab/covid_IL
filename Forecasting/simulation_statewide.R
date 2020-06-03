covid_source('utils.R')

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
  beta_scales,
  frac_underreported,
  seed = NULL,
  format = 'data.frame',
  rprocess_Csnippet,
  initialize = T,
  rinit_Csnippet,
  rmeasure_Csnippet=NULL,
  obsnames = NULL,
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
    input_params, intervention_df, nu_scales, beta_scales, frac_underreported
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
  acc_new_symptomatic <- sprintf("new_symptomatic_infections_%d",c(1:n_age_groups))
  acc_nD <- sprintf("new_deaths_%d",c(1:n_age_groups))
  acc_nHD <- sprintf("new_hosp_deaths_%d",c(1:n_age_groups))
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
    accum_names <- c(accum_names, paste0(acc_new_symptomatic, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_nD, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_nHD, "_", i))
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
  print('Simulating')
  output_df <- simulate(
    nsim = nsim,
    seed = seed,
    times = times,
    t0 = times[1],
    rinit = rinit_Csnippet,
    rprocess = euler(rprocess_Csnippet, delta.t = delta_t),
    rmeasure = rmeasure_Csnippet,
    
    params = params,
    covar = covar_table,
    
    statenames = c(state_names, accum_names),
    paramnames = names(params),
    accumvars = accum_names,
    obsnames = obsnames,
    format = format
  )
  }
  if(initialize == F){
    output_df <- simulate(
      nsim = 1,
      seed = seed,
      times = times,
      t0 = times[1],
      rprocess = euler(rprocess_Csnippet, delta.t = delta_t),
      rmeasure = rmeasure_Csnippet,
      
      params = pfilter_traj_params,
      covar = covar_table,
      
      statenames = c(state_names, accum_names),
      paramnames = names(params),
      accumvars = accum_names,
      obsnames = obsnames,
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
      state = c("E","A","P","IM","IM_dead","IS", "IH1_", "IH2_", "IH3_", "IC2_", "IC3_", "IH4_"),
      index = c(
        alpha_E,
        alpha_A, alpha_P,
        alpha_IM, alpha_IM,
        alpha_IS,
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
      "region_to_test" = region_to_test,
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
      "theta_test"=theta_test,
      "nu_1"= nu_1,
      "nu_3"= nu_3,
      "nu_m"= nu_m,
      "dispersion" = dispersion,
      "frac_severe_init" = frac_severe_init,
      "t_reporting_adjustment" = t_reporting_adjustment,
      "lower_bound_reporting_uncertainty" = lower_bound_reporting_uncertainty,
      "beta_noise_amplitude" = beta_noise_amplitude,
      "b_elderly" = b_elderly
    )
    
    for(i in c(1:length(population_list))){
      params[paste0(paste0("N_",c(1:n_age_groups)),"_",i)] = population_list[[i]]$POPULATION
    }
    
    for(i in c(1:length(population_list))){
      params[paste0(paste0("age_dist_",c(1:n_age_groups)),"_",i)] = population_list[[i]]$POPULATION / sum(population_list[[i]]$POPULATION)
    }
    
    # Beta 2's by region
    if (n_regions == 3){
      params[paste0("beta2_",c(1:n_regions))] = c(beta2_1, beta2_2, beta2_3)
      params[paste0("beta1_",c(1:n_regions))] = c(beta1_1, beta1_2, beta1_3)
    
      # Num init by region 
      params[paste0("num_init_",c(1:n_regions))] = c(num_init_1, num_init_2, num_init_3)
    } else if(n_regions == 4){
      params[paste0("beta2_",c(1:n_regions))] = c(beta2_1, beta2_2, beta2_3, beta2_4)
      params[paste0("beta1_",c(1:n_regions))] = c(beta1_1, beta1_2, beta1_3, beta1_4)
      params[paste0("region_non_hosp_", 1:n_regions)] = c(region_non_hosp_1, region_non_hosp_2, region_non_hosp_3, region_non_hosp_4)
      # Num init by region 
      params[paste0("num_init_",c(1:n_regions))] = c(num_init_1, num_init_2, num_init_3, num_init_4)      
    }
    
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
    params[paste0("mu_m_",c(1:n_age_groups))] = rep(1/inv_mu_m, n_age_groups)
    
    # Probabilities
    params[paste0("rho_",c(1:n_age_groups))] = unlist(
      input_params[paste0('rho_', c(1:n_age_groups))]
    )
    params[paste0("phi_",c(1:n_age_groups))] = unlist(
      input_params[paste0('phi_', c(1:n_age_groups))]
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
    params[paste0("psi4_",c(1:n_age_groups))] = unlist(
      input_params[paste0('psi4_', c(1:n_age_groups))]
    )
    params[paste0("q_",c(1:n_age_groups))] = unlist(
      input_params[paste0('q_', c(1:n_age_groups))]
    )

    params[paste0("age_beta_scales_",c(1:n_age_groups))] = unlist(
      input_params[paste0('age_beta_scales_', c(1:n_age_groups))]
    )

  for(k in c(1:n_regions)){
    home_contact_vec_region <- array()
    for(i in c(1:n_age_groups)){
      for(j in c(1:n_age_groups)){
        home_contact_vec_region[paste0("C_home_",j,"_",i)] = 0
      }
    }
    home_contact_vec_region = home_contact_vec_region[!is.na(home_contact_vec_region)]
    names(home_contact_vec_region) <- paste0(names(home_contact_vec_region), "_",k)
    params = c(params, home_contact_vec_region)
  }
    
    for(k in c(1:n_regions)){
      work_contact_vec_region <- array()
      for(i in c(1:n_age_groups)){
        for(j in c(1:n_age_groups)){
          work_contact_vec_region[paste0("C_work_",j,"_",i)] = 0
        }
      }
      work_contact_vec_region = work_contact_vec_region[!is.na(work_contact_vec_region)]
      names(work_contact_vec_region) <- paste0(names(work_contact_vec_region), "_",k)
      params = c(params, work_contact_vec_region)
    }
    
    for(k in c(1:n_regions)){
      school_contact_vec_region <- array()
      for(i in c(1:n_age_groups)){
        for(j in c(1:n_age_groups)){
          school_contact_vec_region[paste0("C_school_",j,"_",i)] = 0
        }
      }
      school_contact_vec_region = school_contact_vec_region[!is.na(school_contact_vec_region)]
      names(school_contact_vec_region) <- paste0(names(school_contact_vec_region), "_",k)
      params = c(params, school_contact_vec_region)
    }
    
    for(k in c(1:n_regions)){
      other_contact_vec_region <- array()
      for(i in c(1:n_age_groups)){
        for(j in c(1:n_age_groups)){
          other_contact_vec_region[paste0("C_other_",j,"_",i)] = 0
        }
      }
      other_contact_vec_region = other_contact_vec_region[!is.na(other_contact_vec_region)]
      names(other_contact_vec_region) <- paste0(names(other_contact_vec_region), "_",k)
      params = c(params, other_contact_vec_region)
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
    
    intervention_df = data.frame(
        t_start = c(t_start), # t_start and t_end are vectors that have one value for each distinct intervention
        t_end = c(t_end))   
      for(i in c(1:n_regions)){
        intervention_df[,paste0("scale_work_",i)] <- scalework[[i]] ## scalework is a list of length n_regions. Each list item is a region-specific vector with one value for each distinct intervention. 
      }
      for(i in c(1:n_regions)){
        intervention_df[,paste0("scale_home_",i)] <- scalehome[[i]]
      }
      for(i in c(1:n_regions)){
        intervention_df[,paste0("scale_school_",i)] <- scaleschool[[i]]
      }
      for(i in c(1:n_regions)){
        intervention_df[,paste0("scale_other_",i)] <- scaleother[[i]]
      }
    intervention_df
  }) -> intervention_df
  intervention_df
}

simulate_pomp_covid__init_covariate_table <- function(input_params, intervention_df, nu_scales, beta_scales, frac_underreported) {
  n_interventions <- nrow(intervention_df)
  
    # Specify interventions via covariate table 
    name_vec <- names(intervention_df %>% select(contains("scale")))
    
    covar_table_interventions <- data.frame(time = c(1:input_params$tmax),
                                            use_post_intervention_beta = 0)
    for(i in c(1:length(name_vec))){
      covar_table_interventions[,name_vec[i]] = 1
    }
    for(i in c(1:nrow(intervention_df))){
      for(j in c(1:length(name_vec))){
        covar_table_interventions[covar_table_interventions$time >= intervention_df[i,]$t_start & covar_table_interventions$time < intervention_df[i,]$t_end, names(covar_table_interventions) == name_vec[j]] = intervention_df[i, names(intervention_df) == name_vec[j]] 
      }
    }
  covar_table_interventions[covar_table_interventions$time >= intervention_df[1,]$t_start, 'use_post_intervention_beta'] = 1
  
  covar_table_interventions$nu_scale = nu_scales[covar_table_interventions$time, 'nu_scale']
  
  covar_table_interventions$add_noise_to_beta = beta_scales[covar_table_interventions$time, 'add_noise_to_beta']
  
  for (region in 1:input_params$n_regions){
      col = paste0('frac_underreported_',region)
      col_se = paste0('frac_underreported_se_',region)
      covar_table_interventions[[col]] = frac_underreported[covar_table_interventions$time, col]
      covar_table_interventions[[col_se]] = frac_underreported[covar_table_interventions$time, col_se]
  }

  for (region in 1:input_params$n_regions){
      col = paste0('scale_beta_',region)
      covar_table_interventions[[col]] = beta_scales[covar_table_interventions$time, col]
  }
  covar_table_interventions
}

#' Make output table with the following columns:
#' Simulation ID, Time, Compartment, Age group, Cases
 
process_pomp_covid_output <- function(sim_result, agg_regions=T) {
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
  df_sim_output %>% mutate(Region=substr(Compartment, 
    nchar(as.character(Compartment)), 
    nchar(as.character(Compartment))),
    Compartment = case_when(
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
      startsWith(as.character(Compartment), "new_hosp_deaths") ~ "nHD",
      startsWith(as.character(Compartment), "new_mild") ~ "nM",
      startsWith(as.character(Compartment), "new_symptomatic") ~ "nS",
      startsWith(as.character(Compartment), "ObsHospDeaths") ~ "Reported hd",
      startsWith(as.character(Compartment), "ObsNonHospDeaths") ~ "Reported nhd",
      startsWith(as.character(Compartment), "ObsICU") ~ "ObsICU"
  )
  ) -> df_sim_output

  names(df_sim_output) = c('SimID', 'Time', 'Compartment', 'Cases', 'Region')

  if (agg_regions){
      df_sim_output %>%
        group_by(SimID, Time, Compartment) %>%
        summarize(Cases = sum(Cases)) -> df_sim_output
    } else{
      df_sim_output %>%
        group_by(SimID, Time, Compartment, Region) %>%
        summarize(Cases = sum(Cases)) -> df_sim_output
    }

  c(
    sim_result,
    list(plotting_output = df_sim_output)
  )
}

format_for_covid_hub <- function(plotting_output,
  FIPS_code=17,
  forecast_date=Sys.Date()){

  library(purrr)
  library(tidyverse)
  get_target_name <- function(date, forecast_date, compartment){
      days_ahead = as.numeric(date - forecast_date)
      name = sapply(compartment, FUN=function(x){ifelse(x=='D', 'cum death', ifelse(x=='nD', 'inc death', 'inc hosp'))})
      target = sprintf("%s day ahead %s", days_ahead, name)
      return(target)
  }

  p <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  p_names <- map_chr(p, ~.x)
  p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
    set_names(nm = p_names)
  p_funs$`NA` = function(.x){mean(.x)}

  plotting <- as.data.frame(plotting_output)
  plotting <- plotting %>% mutate(Time=as.Date('2020-01-14') + Time)
  quantiles = plotting %>% group_by(Time, Compartment) %>% 
      summarize_at(vars(Cases), funs(!!!p_funs)) %>%
      ungroup() %>%
      filter(Time > forecast_date, Compartment %in% c('D', 'nD', 'H'))

  final_frame <- quantiles %>%     
    gather('quantile', 'value', 3:ncol(quantiles)) %>%
    mutate(location=FIPS_code,
             location_name='Illinois',
             forecast_date=forecast_date,
             type=case_when((quantile=='NA') ~'point',
                            (quantile!='NA') ~'quantile'),
             target=get_target_name(Time, forecast_date, Compartment),
             target_end_date=as.character(Time)
            ) %>%
    select('forecast_date','target','target_end_date','location','location_name','type','quantile','value')
  final_frame
}

format_for_covid_hub_week <- function(plotting_output,
  FIPS_code=17,
  forecast_date=Sys.Date()){

  library(purrr)
  library(tidyverse)
  get_target_name <- function(epi_week, epi_start_week, compartment){
      days_ahead = epi_week - epi_start_week
      name = sapply(compartment, FUN=function(x){ifelse(x=='D', 'cum death', ifelse(x=='nD', 'inc death', 'inc hosp'))})
      target = sprintf("%s wk ahead %s", days_ahead, name)
      return(target)
  }

  p <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  p_names <- map_chr(p, ~.x)
  p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
    set_names(nm = p_names)
  p_funs$`NA` = function(.x){mean(.x)}

  plotting <- as.data.frame(plotting_output)%>% 
      mutate(Time=as.Date('2020-01-14') + Time, 
            Epiweek=MMWRweek(Time)$MMWRweek,
            weekday=weekdays(Time))
  
  incident_deaths <- plotting  %>% 
      filter(Compartment == 'nD') %>% 
      group_by(SimID, parset, Epiweek) %>% summarize(Cases=sum(Cases), Compartment='nD') %>% 
      select(SimID, parset, Epiweek, Compartment, Cases) %>% ungroup()

  cumulative_deaths <-  plotting %>% 
      filter(Compartment =='D', weekday=='Saturday') %>% 
      select(SimID, parset, Epiweek, Compartment, Cases)
  
  plotting = bind_rows(incident_deaths, cumulative_deaths)
  
  
   if (weekdays(forecast_date) %in% c('Sunday','Monday')){
      epi_start_week = MMWRweek(forecast_date)$MMWRweek - 1
  } else{
      epi_start_week = MMWRweek(forecast_date)$MMWRweek
  }
  
  quantiles = plotting %>% group_by(Epiweek, Compartment) %>% 
      summarize_at(vars(Cases), funs(!!!p_funs)) %>%
      ungroup() %>% 
      filter(Epiweek > epi_start_week, Compartment %in% c('D', 'nD'))
  


  final_frame <- quantiles %>%     
    gather('quantile', 'value', 3:ncol(quantiles)) %>%
    mutate(location=FIPS_code,
             location_name='Illinois',
             forecast_date=forecast_date,
             type=case_when((quantile=='NA') ~'point',
                            (quantile!='NA') ~'quantile'),
             target=get_target_name(Epiweek, epi_start_week, Compartment),
             target_end_date=sapply(Epiweek, FUN=function(x){as.character(MMWRweek2Date(2020, x, 7))})
            ) %>%
    select('forecast_date','target','target_end_date','location','location_name','type','quantile','value')
    final_frame
}

add_non_hospitalized_deaths = function(sim_full,
                                       pars,
                                       regional_aggregation=T,
                                       n_age_groups=9,
                                       n_regions=3){
  
  sim_raw_deaths_total <- sim_full$raw_simulation_output %>% select(contains(c("new_deaths")))
  sim_raw_deaths_hosp <- sim_full$raw_simulation_output %>% select(contains(c("new_hosp")))
  
  df_new_nonhosp_deaths <- data.frame(Time =sim_full$raw_simulation_output$time) %>%
    mutate(SimID = sim_full$raw_simulation_output$.id)
  
  # Make sure to loop over each region
  for (region in c(1:n_regions))
  {    
    for(i in c(1:n_age_groups)){
      df_new_nonhosp_deaths[,sprintf("NHD_%s_%s",i, region)] <- sim_raw_deaths_total[,sprintf('new_deaths_%s_%s', i, region)] - sim_raw_deaths_hosp[,sprintf('new_hosp_deaths_%s_%s', i, region)]
    }
  }
  df_new_nonhosp_deaths <- df_new_nonhosp_deaths %>%
    melt(id.vars = c("SimID", "Time")) %>%
    rename(
      Compartment = variable,
      Cases = value
    )
  df_new_nonhosp_deaths <-  df_new_nonhosp_deaths %>% mutate(Region=substr(Compartment, 
                                                                           nchar(as.character(Compartment)), 
                                                                           nchar(as.character(Compartment))),
                                                             Compartment = case_when(
                                                               startsWith(as.character(Compartment), "NHD") ~ "nDnH"))
  
  
  if(regional_aggregation){
    df_new_nonhosp_deaths <- df_new_nonhosp_deaths %>%
      group_by(SimID, Time, Compartment) %>% 
      summarize(Cases=sum(Cases)) %>% 
      ungroup()
    
    df_cum = df_new_nonhosp_deaths %>%
      group_by(SimID) %>% 
      arrange(Time) %>% 
      mutate(Cases = cumsum(Cases), 
             Compartment='DnH') %>% ungroup()
    df_non_hosp = bind_rows(df_new_nonhosp_deaths, df_cum)
    
  } else{
    
    df_cum = df_new_nonhosp_deaths %>%
      group_by(SimID, Region) %>% 
      arrange(Time) %>% 
      mutate(Cases = cumsum(Cases),
             Compartment='DnH') %>% ungroup()
    df_non_hosp = bind_rows(df_new_nonhosp_deaths, df_cum)
  }
  return(df_non_hosp)
}

get_reported_non_hospitalized_deaths = function(df_input,
                                                regional_aggregation,
                                                lower_bound_reporting = 0.25,
                                                upper_bound_reporting = 0.75){
  if(regional_aggregation){
     df_nHD <- df_input %>% 
    filter(Compartment == "nDnH") %>% 
    group_by(parset, SimID) %>% 
    mutate(Cases_reported = round(runif(1,lower_bound_reporting,upper_bound_reporting)*Cases)) %>% 
    ungroup() 
    } else{
         df_nHD <- df_input %>% 
          filter(Compartment == "nDnH") %>% 
          group_by(parset, SimID, Region) %>% 
    mutate(Cases_reported = round(runif(1,lower_bound_reporting,upper_bound_reporting)*Cases)) %>% 
    ungroup()   
    }
  return(df_nHD)
}

get_scale = function(t_logistic_start,
                     intervention_lift,
                     simstart,
                     simend,
                     max_scales){
    
    times = seq(1, simend, 1)
    
    logistic = function(x, mscale, shift=intervention_lift){
        mean = (mscale-1)/(1+exp(-(x - shift))) + 1
        mean
    }

    raw_scales = data.frame(time=times, 
                            scale_beta_1=logistic(times, mscale=max_scales[1]),
                            scale_beta_2=logistic(times, mscale=max_scales[2]),
                            scale_beta_3=logistic(times, mscale=max_scales[3]),
                            scale_beta_4=logistic(times, mscale=max_scales[4])) %>%
        mutate(scale_beta_1 = case_when((time < t_logistic_start) ~1,
                                 (time>=t_logistic_start) ~scale_beta_1),
                scale_beta_2 = case_when((time < t_logistic_start) ~1,
                                                 (time>=t_logistic_start) ~scale_beta_2),
                scale_beta_3 = case_when((time < t_logistic_start) ~1,
                                                 (time>=t_logistic_start) ~scale_beta_3),
                scale_beta_4 = case_when((time < t_logistic_start) ~1,
                                                 (time>=t_logistic_start) ~scale_beta_4),
               add_noise_to_beta = case_when((time < t_logistic_start) ~ 0,
                                             (time >= t_logistic_start) ~ 1))
    
    row.names(raw_scales) = raw_scales$time
    raw_scales
    
}

add_interventions = function(intervention_file, parameters){
  stopifnot(require(foreach))
  pars <- parameters

  # Add the intervention starts and ends
  interventions = read.csv(intervention_file) %>% arrange(t_start, region)
  pars$t_start = unique(interventions$t_start)
  pars$t_end = unique(interventions$t_end)
  
  # Add in scalings
  interventions = interventions %>% group_by(region, t_start) %>% 
    summarise(scalework=median(scalework), scaleschool=median(scaleschool), scalehome=median(scalehome), scaleother=median(scaleother)) %>%
    ungroup()
  for(location in c('scalework','scaleschool','scalehome','scaleother')){
    foreach(r=1:pars$n_regions) %do% {
    temp = interventions %>% filter(region==r) %>% arrange(t_start)
    temp[[location]]
    } -> scales
    pars[[location]] = scales
  }
  return(pars)
}


