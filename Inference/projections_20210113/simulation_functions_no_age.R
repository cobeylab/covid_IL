covid_source('utils.R')
library(stringr)

#' Simulate COVID pomp model
simulate_pomp_covid <- function(
  n_regions,
  input_params, 
  delta_t, 
  population_list, 
  death_reporting_covar, # covariates to calculate death reporting
  nonhosp_death_covar, # covariates to calculate fraction of deaths non-hospitalized
  emr_report_covar,
  beta_covar=NULL, # covariate of raw beta values
  beta_scale_covar=NULL, # covariate of values to scale beta, e.g., mobility
  seed = NULL,
  format = 'data.frame',
  time_column='time',
  rprocess_Csnippet,
  initialize = T,
  inference=F,
  rinit_Csnippet,
  n_age_groups=1,
  nsim=NULL, 
  rmeasure_Csnippet=NULL,
  dmeasure_Csnippet=NULL,
  obsnames = NULL,
  initial_params = NULL,
  data=NULL,
  transformations=NULL
) {
  library(pomp)
  
  # Initialize states
  subcompartment_df <- simulate_pomp_covid__init_subcompartment_df(input_params)

  # Set up parameters for pomp
  params <- simulate_pomp_covid__init_parameters(input_params, population_list)

  # Set up covariate table for pomp
  covar_table_interventions <- simulate_pomp_covid__init_covariate_table(
    input_params, 
    beta_scale_covar, 
    beta_covar, 
    death_reporting_covar,
    nonhosp_death_covar, 
    emr_report_covar
  )

  covar_table <- covariate_table(
    covar_table_interventions,
    order = "constant",
    times = "time"
  )

  # Actually call pomp
  state_names <- simulate_pomp_covid__init_state_names(n_age_groups, n_regions, subcompartment_df)
  ## Set up accumulator variables by age group and then by regions
  acc_nD <- sprintf("new_deaths_%d",c(1:n_age_groups))
  acc_nHD <- sprintf("new_hosp_deaths_%d",c(1:n_age_groups))
  acc_nNHD <- sprintf("new_nonhosp_deaths_%d", c(1:n_age_groups))
  acc_nH <- sprintf("new_hospitalizations_%d", c(1:n_age_groups))
  acc_INC <- sprintf("Inc_%d", c(1:n_age_groups))
  
  accum_names <- array()
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_nD, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_nHD, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_nNHD, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_nH, "_", i))
  }
  for(i in c(1:n_regions)){
    accum_names <- c(accum_names, paste0(acc_INC, "_", i))
  }
  
  accum_names <- accum_names[!is.na(accum_names)]
  times <- c(ceiling(input_params$tmin):ceiling(input_params$tmax)) 
  if(initialize == T & inference==F){
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
      format = format)

    return_value = list(
      n_age_groups = n_age_groups,
      raw_simulation_output = output_df,
      params = params,
      state_names = state_names,
      interventions = covar_table_interventions
    )
  } else if(inference == T){
      return_value =  pomp(
      times = time_column,
      t0 = times[1],
      data=data,
      partrans=transformations,
      rinit = rinit_Csnippet,
      rprocess = euler(rprocess_Csnippet, delta.t = delta_t),
      dmeasure = dmeasure_Csnippet,
      params = params,
      covar = covar_table,
      statenames = c(state_names, accum_names),
      paramnames = names(params),
      accumvars = accum_names,
      obsnames = obsnames
      )
    }
  return_value
} 

simulate_pomp_covid__init_subcompartment_df <- function(input_params) {
  with(input_params, {
    data.frame(
      state = c("E","P","IM","IM_dead","IS", "IH1_", "IH2_", "IC_", "preIC_"),
      index = c(
        alpha_E,
        alpha_P,
        alpha_IM, alpha_IM,
        alpha_IS,
        alpha_IH1, alpha_IH2,
        alpha_IC,
        alpha_preIC
      )
    )
  })
}

simulate_pomp_covid__init_state_names <- function(n_age_groups, n_regions, subcompartment_df) {
  ## State variables with variable classes per age group
  state_names_S_age = c(sprintf("S_%d",c(1:n_age_groups)))
  state_names_R_age = c(sprintf("R_%d",c(1:n_age_groups)))
  state_names_D_age <- c(sprintf("D_%d",c(1:n_age_groups)))
  
  state_names_trans_age <- c(sprintf("transmissionRate_%d",c(1:n_age_groups)))
  state_names_transMean_age <- c(sprintf("transmissionMean_%d",c(1:n_age_groups)))
  
  state_names_hfr_age <- c(sprintf("HFRtrack_%d",c(1:n_age_groups)))
  state_names_hfrMean_age <- c(sprintf("HFRMean_%d",c(1:n_age_groups)))

  state_names_icuTrack_age <- c(sprintf("ICUtrack_%d",c(1:n_age_groups)))
  state_names_icuMean_age <- c(sprintf("ICUMean_%d",c(1:n_age_groups)))


  state_names_phi_age <- c(sprintf("phitrack_%d",c(1:n_age_groups)))
  state_names_ihr_age <- c(sprintf("IHRtrack_%d",c(1:n_age_groups)))
  state_names_ifr_age <- c(sprintf("IFRtrack_%d",c(1:n_age_groups)))
  state_names_deathreport_age <- c(sprintf("DeathReportTrack_%d",c(1:n_age_groups)))

  state_names_S <- array()
  state_names_R <- array()
  state_names_D <- array()
  state_names_trans <- array()
  state_names_transMean <- array()
  state_names_hfr <- array()
  state_names_hfrMean <- array()
  state_names_icuTrack <- array()
  state_names_icuMean <- array()
  state_names_phi <- array()
  state_names_ihr <- array()
  state_names_ifr <- array()
  state_names_deathreport <- array()

  for(i in c(1:n_regions)){
    state_names_S <- c(state_names_S, paste0(state_names_S_age, "_", i))
    state_names_R <- c(state_names_R, paste0(state_names_R_age, "_", i))
    state_names_D <- c(state_names_D, paste0(state_names_D_age, "_", i))
    state_names_trans <- c(state_names_trans, paste0(state_names_trans_age, "_", i))
    state_names_transMean <- c(state_names_transMean, paste0(state_names_transMean_age, "_", i))
    
    state_names_hfr <- c(state_names_hfr, paste0(state_names_hfr_age, "_", i))
    state_names_hfrMean <- c(state_names_hfrMean, paste0(state_names_hfrMean_age, "_", i))
    state_names_icuTrack <- c(state_names_icuTrack, paste0(state_names_icuTrack_age, "_", i))
    state_names_icuMean <- c(state_names_icuMean, paste0(state_names_icuMean_age, "_", i))

    state_names_phi <- c(state_names_phi, paste0(state_names_phi_age, "_", i))
    state_names_ihr <- c(state_names_ihr, paste0(state_names_ihr_age, "_", i))
    state_names_ifr <- c(state_names_ifr, paste0(state_names_ifr_age, "_", i))
    state_names_deathreport = c(state_names_deathreport, paste0(state_names_deathreport_age, '_', i))
  }
    
  state_names <- c(state_names_S, state_names_R, state_names_D, 
    state_names_trans, state_names_transMean,
    state_names_hfr, state_names_hfrMean,
    state_names_icuTrack, state_names_icuMean,
    state_names_phi, state_names_ihr, state_names_ifr,
    state_names_deathreport
    )

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
  input_params, population_list
) {
  with(input_params, {
    # Set up parameters for model
    params = c(
      # Control params
      "n_regions" = n_regions,
      "region_to_test" = region_to_test,
      "eta"=eta,
      "region_specific_HFR"=region_specific_HFR,
      "alpha_E"= alpha_E,
      "alpha_P"= alpha_P,
      "alpha_IM"= alpha_IM,
      "alpha_IS"= alpha_IS,
      "alpha_IH1" = alpha_IH1, 
      "alpha_IH2" = alpha_IH2,
      "alpha_IC" = alpha_IC,
      "alpha_preIC" = alpha_preIC,      
      # Constant rates
      "sigma" = 1/inv_sigma,
      "zeta_s" = 1/inv_zeta_s,
      "mu_m" = 1/inv_mu_m,
      "gamma_m" = 1/inv_gamma_m,
      "zeta_h" = 1/inv_zeta_h,
      "pre_icu_rate"=1/inv_pre_icu_rate,

      # Time-varying params
      "mu_min"=mu_min,
      "mu_max"=mu_max,
      "ICU_min"=ICU_min,
      "ICU_max"=ICU_max,
      "gamma_min"=gamma_min,
      "gamma_max"=gamma_max,
      "zeta_icu_min"=zeta_icu_min,
      "zeta_icu_max"=zeta_icu_max,
      "HFR_min"=HFR_min,
      "HFR_max"=HFR_max,
      "hosp_t_min"=hosp_t_min,
      "hosp_t_max"=hosp_t_max,
      "IFR_constraint"=IFR_constraint
    )
    
    # Region-specific parameters
    params[paste0("num_init_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('num_init_', c(1:n_regions))])
      )
    params[paste0("n_changepoints_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('n_changepoints_', c(1:n_regions))])
      )
    params[paste0("n_HFR_changepoints_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('n_HFR_changepoints_', c(1:n_regions))])
      )
    params[paste0("n_ICU_changepoints_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('n_ICU_changepoints_', c(1:n_regions))])
      )
    params[paste0("phi_scale_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('phi_scale_', c(1:n_regions))])
      )
    params[paste0("inv_mu_0_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('inv_mu_0_', c(1:n_regions))])
      )
    params[paste0("inv_gamma_0_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('inv_gamma_0_', c(1:n_regions))])
      )
    params[paste0("inv_zeta_icu_0_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('inv_zeta_icu_0_', c(1:n_regions))])
      )
    params[paste0("inv_zeta_icu_f_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('inv_zeta_icu_f_', c(1:n_regions))])
      )
    params[paste0("inv_gamma_f_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('inv_gamma_f_', c(1:n_regions))])
      )
    params[paste0("inv_mu_f_",c(1:n_regions))] = unlist(
        unlist(input_params[paste0('inv_mu_f_', c(1:n_regions))])
      )

    # Region and changepoint specific
    for (region in 1:n_regions){
      eval_str = paste0('n_changepoints_',region)
      n_changepoints_reg = eval(parse(text=eval_str))
      changepoint_names = paste0("changepoints_",c(1:n_changepoints_reg),'_',region)
      params[changepoint_names] = unlist(unlist(input_params[changepoint_names]))

    }
    for (region in 1:n_regions){
      eval_str = paste0('n_changepoints_',region)
      n_changepoints_reg = eval(parse(text=eval_str))
      beta_names = paste0("beta_values_",c(1:n_changepoints_reg),'_',region)
      params[beta_names] = unlist(unlist(input_params[beta_names]))
    }

    for (region in 1:n_regions){
      eval_str = paste0('n_HFR_changepoints_',region)
      n_changepoints_reg = eval(parse(text=eval_str))
      changepoint_names = paste0("HFR_changepoint_values_",c(1:n_changepoints_reg),'_',region)
      params[changepoint_names] = unlist(unlist(input_params[changepoint_names]))
    }
    for (region in 1:n_regions){
      eval_str = paste0('n_HFR_changepoints_',region)
      n_changepoints_reg = eval(parse(text=eval_str))
      hfr_names = paste0("HFR_values_",c(1:n_changepoints_reg),'_',region)
      params[hfr_names] = unlist(unlist(input_params[hfr_names]))
    }


    for (region in 1:n_regions){
      eval_str = paste0('n_ICU_changepoints_',region)
      n_changepoints_reg = eval(parse(text=eval_str))
      changepoint_names = paste0("ICU_changepoint_values_",c(1:n_changepoints_reg),'_',region)
      params[changepoint_names] = unlist(unlist(input_params[changepoint_names]))
    }
    for (region in 1:n_regions){
      eval_str = paste0('n_ICU_changepoints_',region)
      n_changepoints_reg = eval(parse(text=eval_str))
      hfr_names = paste0("ICU_values_",c(1:n_changepoints_reg),'_',region)
      params[hfr_names] = unlist(unlist(input_params[hfr_names]))
    }

    params[paste0("N_",c(1:n_regions))] = population_list %>%
      group_by(covid_region) %>%
      dplyr::summarize(pop = sum(POPULATION)) %>%
      ungroup() %>%
      arrange(as.numeric(covid_region)) %>%
      select(pop) %>%
      unlist(use.names=F)


    params
  })
}

simulate_pomp_covid__init_covariate_table <- function(input_params, 
  beta_scale_covar, 
  beta_covar, 
  death_reporting_covar,
  nonhosp_death_covar, 
  emr_report_covar) {
    
  covar_table_interventions <- data.frame(time = c(1:input_params$tmax))
  
  # add death underreporting
  covar_table_interventions = left_join(covar_table_interventions, death_reporting_covar, by='time')

  # add in non-hospitalized death info
  covar_table_interventions = left_join(covar_table_interventions, nonhosp_death_covar, by='time')

  # add emr_reporting
  covar_table_interventions = left_join(covar_table_interventions, emr_report_covar, by='time')
  covar_table_interventions
}
 
process_pomp_covid_output <- function(sim_result, agg_regions=T) {

  num_age_classes = 1
  num_regions = 11
  max_subcompartments = 3
  # compartments with no subcompartments
  no_sub = c("S", "D", "R", "transmissionRate", "HFRtrack","phitrack","IHRtrack", "IFRtrack", "Inc", "new_deaths", "new_hosp_deaths", "new_nonhosp_deaths", "new_hospitalizations",
    "transmissionMean", "ICUtrack", "ICUMean", "HFRMean", "DeathReportTrack",
   paste0("E", 1:max_subcompartments), 
    paste0("P", 1:max_subcompartments), 
    paste0("IS", 1:max_subcompartments), 
    paste0("IM", 1:max_subcompartments),
    paste0("IM_dead", 1:max_subcompartments)
    )
  subcompartments = c("IH1", "IH2", "IC", "preIC")
  obs_prefix = c("ObsDeaths", "ObsHosp", "ObsHospDeaths","ObsICU")

  nosub_suffixes = expand.grid(no_sub, 1:num_age_classes, 1:num_regions)
  nosub_suffixes$state_name = paste(nosub_suffixes$Var1, nosub_suffixes$Var2, nosub_suffixes$Var3,  sep='_')
  nosub_suffixes$Var1 = gsub('[[:digit:]]+', '', nosub_suffixes$Var1)

  sub_suffixes = expand.grid(subcompartments, 1:max_subcompartments, 1:num_age_classes, 1:num_regions)
  sub_suffixes$state_name = paste(sub_suffixes$Var1, sub_suffixes$Var2, sub_suffixes$Var3, sub_suffixes$Var4,  sep='_')
  sub_suffixes$Var1 = gsub('[[:digit:]]+', '', sub_suffixes$Var1)

  obs_suffixes = expand.grid(obs_prefix, 1:num_regions)
  obs_suffixes$state_name = paste(obs_suffixes$Var1, obs_suffixes$Var2,  sep='_')


  name_map = c(as.vector(nosub_suffixes$Var1), as.vector(sub_suffixes$Var1), as.vector(obs_suffixes$Var1))
  names(name_map) = c(nosub_suffixes$state_name, sub_suffixes$state_name, obs_suffixes$state_name)
  region_map = c(as.vector(nosub_suffixes$Var3), as.vector(sub_suffixes$Var4), as.vector(obs_suffixes$Var2))
  names(region_map) = names(name_map)

  age_map = c(as.vector(nosub_suffixes$Var2), as.vector(sub_suffixes$Var3), as.vector(rep(-1, length(as.vector(obs_suffixes$Var2)))))
  names(age_map) = names(name_map)

  df_sim <- sim_result
  #params <- sim_result$params

  print('Melting')
  df_sim_output = df_sim %>% 
      melt(id.vars = c(".id", "time")) %>% 
      rename(
          Compartment = variable,
          Cases = value
      ) %>% 
      drop_na() %>%
      mutate(
        Compartment=as.character(Compartment),
        Region= region_map[Compartment],
        Compartment = name_map[Compartment]
          ) %>%
      rename(SimID = .id,
        Time = time) %>%
      group_by(Compartment, Region, Time, SimID) %>%
      summarize(Cases = sum(Cases)) %>%
      ungroup() %>%
      # reset initial incidence
      group_by(Region, SimID) %>%
      mutate(Cases = if_else(Compartment == "Inc" & Time == min(Time), NA_real_, Cases)) %>%
      ungroup()
    print('Calculating prevalence')
  prevalenceComps = c('E', 'P', 'IS', 'IM')

  prevdf = df_sim_output %>%
    filter(Compartment %in% prevalenceComps) %>%
    group_by(Region, Time, SimID) %>%
    summarize(Cases=sum(Cases)) %>%
    mutate(Compartment="Total infectious") %>%
    select(names(df_sim_output))

  df_sim_output = bind_rows(df_sim_output, prevdf)

  print("Finished processing")

  list(plotting_output = df_sim_output)

}

run_pfilter_and_output_result <- function(n_reps_pfilter, 
  n_particles_pfilter, 
  jobid,
  arrayid,
  output_dir,
  pomp_object,
  t_mif,
  return_df =F){
  require(foreach)
  chainID = sprintf("%s_%s", jobid, arrayid)
  output_filename = sprintf("%s/output_pfilter_%s_%s.csv", output_dir, jobid, arrayid)
  t = Sys.time()

  ll <- logmeanexp(replicate(n_reps_pfilter,logLik(pfilter(pomp_object,Np=n_particles_pfilter))),se=TRUE)
  runtime = Sys.time() - t

  finaloutput = (data.frame(as.list(pomp_object@params),
                      loglik=ll[1],
                      loglik_se=ll[2],
                      chain = chainID,
                      t_pfilter=runtime,
                      t_mif=t_mif,
                      nparticles=n_particles_pfilter,
                      n_reps=n_reps_pfilter))    

  write.table(finaloutput, file = output_filename, sep = ",",col.names = TRUE, row.names=FALSE, append=TRUE)
  if(return_df == T){
    return(finaloutput)
  }
}

run_mif_and_output_result <- function(
  jobid,
  arrayid,
  output_dir,
  pomp_object,
  n_mif,
  rw_vec,
  cooling_rate,
  n_particles_mif){

  t = Sys.time()

  pomp_object %>% mif2(
    Nmif= n_mif ,
    rw.sd = rw_vec,
    cooling.type="geometric",
    cooling.fraction.50= cooling_rate,
    Np= n_particles_mif
  ) -> mf

  chainID = sprintf("%s_%s", jobid, arrayid)
  output_filename = sprintf("%s/output_mif_%s_%s.rds", output_dir, jobid, arrayid)
  saveRDS(mf, output_filename)
  mf
}

convert_date_to_time = function(d1, dref){
  as.numeric(as.Date(d1) - as.Date(dref))
}


unlogit = function( logit_value, pmax, pmin) {
    (logit_value * (pmax - pmin)) + pmin
}

get_par_value = function( p0, pf, t0, tf, t_now ) {
   slope = (pf - p0) / (tf - t0)
   if ( t_now >= tf){
       value = pf
   } else if (t_now <= t0){
       value = p0
   } else{
       value = p0 + slope * (t_now - t0)
   }
  value
}