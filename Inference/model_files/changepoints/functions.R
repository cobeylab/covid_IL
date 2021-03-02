covid_source('utils.R')
library(stringr)

#' Simulate COVID pomp model
simulate_pomp_covid <- function(
  n_regions,
  input_params, 
  delta_t, 
  population_list, 
  region_covars,
  shared_covars,
  seed = NULL,
  format = 'data.frame',
  time_column='time',
  rprocess_Csnippet,
  initialize = T,
  inference=F,
  rinit_Csnippet,
  global_Csnippet,
  n_age_groups=1,
  nsim=NULL, 
  rmeasure_Csnippet=NULL,
  dmeasure_Csnippet=NULL,
  obsnames = NULL,
  initial_params = NULL,
  data=NULL,
  transformations=NULL,
  use.mle.params = F,
  mle.params=NULL
) {
  # Initialize states
  subcompartment_df <- simulate_pomp_covid__init_subcompartment_df(input_params)

  # Set up parameters for pomp
  if(!use.mle.params){
        params = simulate_pomp_covid__init_parameters(input_params, population_list)
        if (input_params$region_to_test > 0){
            pars_to_keep = c(
              "n_regions",
              "region_to_test",
              "eta",
              "alpha_E",
              "alpha_P",
              "alpha_IM",
              "alpha_IS",
              "alpha_IH1",
              "alpha_IH2",
              "alpha_IC",
              "alpha_preIC",
              "sigma",
              "zeta_s",
              "mu_m",
              "gamma_m",
              "zeta_h",
              "pre_icu_rate",
              "ICU_min",
              "ICU_max",
              "HFR_min",
              "HFR_max",
              "hosp_t_min",
              "hosp_t_max",
              "IFR_constraint",
              "t_data_max",
              "reporting_cap"
            )

            reg_pars = params[grep(sprintf(".*_%s$", input_params$region_to_test), names(params))]
            names(reg_pars) = sub(sprintf("_%s$", input_params$region_to_test), "_1", names(reg_pars))
            params = c(params[pars_to_keep], reg_pars)
            n_regions = 1
            params[['n_regions']] = 1
          }
    } else{
        params = mle.params
        n_regions = 1
        params[['n_regions']] = 1
    }

  # Set up covariate table for pomp
  covar_table_interventions <- simulate_pomp_covid__init_covariate_table(
    input_params, 
    region_covars,
    shared_covars
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
  for(i in c(1:params[['n_regions']] )){
    accum_names <- c(accum_names, paste0(acc_nD, "_", i))
  }
  for(i in c(1:params[['n_regions']] )){
    accum_names <- c(accum_names, paste0(acc_nHD, "_", i))
  }
  for(i in c(1:params[['n_regions']] )){
    accum_names <- c(accum_names, paste0(acc_nNHD, "_", i))
  }
  for(i in c(1:params[['n_regions']] )){
    accum_names <- c(accum_names, paste0(acc_nH, "_", i))
  }
  for(i in c(1:params[['n_regions']] )){
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
      format = format,
      globals = global_Csnippet)

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
      obsnames = obsnames,
      globals = global_Csnippet
      )
    }
  return_value
} 

simulate_pomp_covid__init_subcompartment_df <- function(input_params) {
  with(input_params, {
    data.frame(
      state = c("E","P","IM","IM_dead","IS", "IH1_", "IH2_", "IC_", "preIC_","Evar","Pvar","IMvar","IM_deadvar","ISvar"),
      index = c(
        alpha_E,
        alpha_P,
        alpha_IM, alpha_IM,
        alpha_IS,
        alpha_IH1, alpha_IH2,
        alpha_IC,
        alpha_preIC,
        alpha_E,
        alpha_P,
        alpha_IM, alpha_IM,
        alpha_IS
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
  state_names_hfr_age <- c(sprintf("HFRtrack_%d",c(1:n_age_groups)))
  state_names_icuTrack_age <- c(sprintf("ICUtrack_%d",c(1:n_age_groups)))
  state_names_ihr_age <- c(sprintf("IHRtrack_%d",c(1:n_age_groups)))
  state_names_ifr_age <- c(sprintf("IFRtrack_%d",c(1:n_age_groups)))
  state_names_deathreport_age <- c(sprintf("DeathReportTrack_%d",c(1:n_age_groups)))
  state_names_freqtrack_age <- c(sprintf("freqtrack_%d",c(1:n_age_groups)))

  state_names_S <- array()
  state_names_R <- array()
  state_names_D <- array()
  state_names_trans <- array()
  state_names_hfr <- array()
  state_names_icuTrack <- array()
  state_names_ihr <- array()
  state_names_ifr <- array()
  state_names_deathreport <- array()
  state_names_freqtrack <- array()

  for(i in c(1:n_regions)){
    state_names_S <- c(state_names_S, paste0(state_names_S_age, "_", i))
    state_names_R <- c(state_names_R, paste0(state_names_R_age, "_", i))
    state_names_D <- c(state_names_D, paste0(state_names_D_age, "_", i))
    state_names_trans <- c(state_names_trans, paste0(state_names_trans_age, "_", i))
    
    state_names_hfr <- c(state_names_hfr, paste0(state_names_hfr_age, "_", i))
    state_names_icuTrack <- c(state_names_icuTrack, paste0(state_names_icuTrack_age, "_", i))

    state_names_ihr <- c(state_names_ihr, paste0(state_names_ihr_age, "_", i))
    state_names_ifr <- c(state_names_ifr, paste0(state_names_ifr_age, "_", i))
    state_names_deathreport = c(state_names_deathreport, paste0(state_names_deathreport_age, '_', i))
    state_names_freqtrack = c(state_names_freqtrack, paste0(state_names_freqtrack_age, '_', i))
  }
    
  state_names <- c(state_names_S, state_names_R, state_names_D, 
    state_names_trans,
    state_names_hfr,
    state_names_icuTrack,
    state_names_ihr, state_names_ifr,
    state_names_deathreport,
    state_names_freqtrack
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
      "ICU_min"=ICU_min,
      "ICU_max"=ICU_max,
      "HFR_min"=HFR_min,
      "HFR_max"=HFR_max,
      "hosp_t_min"=hosp_t_min,
      "hosp_t_max"=hosp_t_max,
      "IFR_constraint"=IFR_constraint,
      "t_data_max"=t_data_max,
      "reporting_cap"=reporting_cap
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
  region_covars,
  shared_covars) {
    
  covar_table_interventions <- data.frame(time = c(1:500))
  # add  region-specific covars
  covar_table_interventions = left_join(covar_table_interventions, region_covars, by='time')
  if(input_params$region_to_test > 0){
    times = covar_table_interventions$time
    cols = names(covar_table_interventions)
    final_cols = cols[grep(sprintf(".*_%s$", input_params$region_to_test), cols)]
    covar_table_interventions = covar_table_interventions %>%
      select(all_of(final_cols))
    names(covar_table_interventions) = sub(sprintf("_%s$", input_params$region_to_test), "_1", names(covar_table_interventions))
    covar_table_interventions = covar_table_interventions %>% mutate(time=times)
  }

  # add non-region-specific covars
  covar_table_interventions = left_join(covar_table_interventions, shared_covars, by='time')
  covar_table_interventions
}
 
 
process_pomp_covid_output <- function(sim_result, agg_regions=T) {

  num_age_classes = 1
  num_regions = 11
  max_subcompartments = 3
  # compartments with no subcompartments
  no_sub = c("S", "D", "R", "transmissionRate", "HFRtrack","IHRtrack", "IFRtrack", "Inc", "new_deaths", "new_hosp_deaths", "new_nonhosp_deaths", "new_hospitalizations",
    "ICUtrack", "DeathReportTrack", "freqtrack",
   paste0("E", 1:max_subcompartments), 
    paste0("P", 1:max_subcompartments), 
    paste0("IS", 1:max_subcompartments), 
    paste0("IM", 1:max_subcompartments),
    paste0("IM_dead", 1:max_subcompartments),
   paste0("Evar", 1:max_subcompartments), 
    paste0("Pvar", 1:max_subcompartments), 
    paste0("ISvar", 1:max_subcompartments), 
    paste0("IMvar", 1:max_subcompartments),
    paste0("IM_deadvar", 1:max_subcompartments)
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
  prevalenceComps = c('E', 'P', 'IS', 'IM', 'IM_dead', 'Evar', 'Pvar', 'ISvar', 'IMvar', 'IM_deadvar', 'ISvar')

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

get_rw_sd_str_time_var = function(time_pars, par_names, sd.val, maxt=1000){

  foreach(i=1:(length(time_pars))) %do%{
      t1 = ifelse(i == 1, 0, time_pars[i - 1])
      t2 = ifelse(i == length(time_pars), maxt, time_pars[i + 1])
      parname = par_names[i]
      s = sprintf('%s=ifelse(time >= %s & time <= %s, %s, 0)', parname, t1, t2, sd.val)
      s
  } -> strings
  paste0(as.character(unlist(strings)), sep=', ', collapse='')
}



get_observed_ifr = function(df.in){
  D = df.in %>% ungroup() %>%  
      filter(Compartment=="D") %>%
      group_by(SimID, Time) %>%
      summarize(Dead=sum(Cases)) %>%
      ungroup() %>%
      mutate(Time = Time - 19)
  R = df.in %>% ungroup() %>%  
      filter(Compartment=="R") %>%
      group_by(SimID, Time) %>%
      summarize(Recovered=sum(Cases)) %>%
      ungroup() %>%
      mutate(Time = Time - 9)
  full_join(D, R, by=c('SimID', 'Time')) %>%
      drop_na() %>%
      mutate(IFR = Dead /(Dead+Recovered)) %>%
      select(SimID, Time, IFR) %>%
      gather(Compartment, Cases, IFR) %>%
      mutate(Date=as.Date('2020-01-14') + Time, 
      Source='Simulation') %>%
      drop_na() %>%
      filter(Date>= as.Date('2020-03-20'))
}

get_cumulative = function(df.in, compartment.name, new.compartment.name, source.name){
  df.in %>%  
    filter(Compartment == compartment.name) %>%
    group_by(SimID) %>%
    arrange(Date) %>%
    mutate(Cases=cumsum(Cases)) %>%
    ungroup() %>%
    mutate(Compartment = new.compartment.name) %>%
    mutate(Source=source.name)
}

get_Rt = function(df.in, new_pars, region_to_test){
    source('R0_functions.R')
    with(as.list(new_pars), {
      beta_init = unlogit(beta_values_1_1, 0.9, 0.4)
      HFR_init = unlogit(HFR_values_1_1, HFR_max, HFR_min)
      phi = unlogit(phi_scale_1, IFR_constraint, 0)
      IHR_init = (IFR_constraint - phi) / (HFR_init - phi)
      paras = list(sigma = 1/3.5,
                   zeta_s = 1/1.5,
                   mu_m = 1/15,
                   gamma_m = 1/3.5,
                   zeta_h = 1/5.6,
                   phi = phi)
      r0pop = read.csv(covid_get_path('Data/covid_region_populations.csv')) %>% 
          summarize(POPULATION=sum(POPULATION))
      R0_val = get_R0(region_cons=region_to_test, 
             beta_value=beta_init,
             pop=r0pop,
             paras=paras,
             t=47,
             IHR =IHR_init)

      Rt_scaling_constant = R0_val / beta_init
      list(rtdf = df.in %>%
        filter(Compartment == 'Transmission rate' | Compartment == 'Fraction recovered') %>%
        spread(Compartment, Cases) %>%
        mutate(Rt = Rt_scaling_constant * `Transmission rate` * (1 - `Fraction recovered`)) %>%
        select(-`Transmission rate`, -`Fraction recovered`) %>%
        rename(Cases=Rt) %>%
        mutate(Compartment='R(t)'),
        R0=R0_val)

  })
}

make_base_plotdf = function(plotting_output,
  compartments_to_plot,
  new_compartment_name_map){

  plotting_output %>% ungroup() %>%    
      mutate(Date=as.Date('2020-01-14') + Time) %>%
      filter(Compartment %in% compartments_to_plot) %>%
      mutate(
        Source = ifelse(Compartment %in% c('ObsHospDeaths', 'ObsDeaths', 'ObsHosp', "ObsICU"), 'Simulation (reported)', 'Simulation'),
        Compartment = new_compartment_name_map[Compartment]
        ) %>%
      group_by(Region, Date, Compartment, SimID, Source) %>%
      dplyr::summarize(Cases=sum(Cases)) %>%
      ungroup()
}

plot_generic = function(alloutputs,
  plotting_data,
  new_pars,
  compartments_to_plot,
  new_compartment_name_map,
  plot_order,
  compartment_order,
  loglik,
  popfinal){


  plotout_noage = make_base_plotdf(alloutputs$plotting_output, compartments_to_plot, new_compartment_name_map)
  cumulative_hosp = get_cumulative(plotout_noage, "Reported hospital deaths", 'Reported cumulative hospital deaths', 'Simulation (reported)')
  cumulative_all_deaths = get_cumulative(plotout_noage, "Reported non-hospital deaths", 'Reported cumulative non-hosp deaths', 'Simulation (reported)')
  IFRdf = get_observed_ifr(alloutputs$plotting_output)
  plotout_noage = bind_rows(plotout_noage, IFRdf) %>%
    bind_rows(cumulative_hosp) %>%
    bind_rows(cumulative_all_deaths)
  frac_recovered = plotout_noage %>%
    filter(Compartment == 'Fraction recovered') %>%
    mutate(Cases = Cases / popfinal[region_to_test])
  plotout_noage = bind_rows(plotout_noage %>% filter(Compartment !='Fraction recovered'), frac_recovered)
  Rtoutputs = get_Rt(plotout_noage, new_pars, region_to_test)
  Rt = Rtoutputs$rtdf
  
  rt_model = Rt %>%
    select(-Time) %>%
    drop_na() %>%
    group_by(Date) %>%
    summarize(rt_model_med=quantile(Cases, 0.5, type=3),
              rt_model_lower=quantile(Cases, 0.025, type=3),
              rt_model_upper=quantile(Cases, 0.975, type=3)) %>%
    mutate(geography_modeled = sprintf('covidregion_%s', region_to_test),
           scenario_name = 'baseline') %>%
    ungroup() %>%
    rename(date=Date)

  write.csv(rt_model, sprintf('%s/model_rt_%s_%s.csv', human_output_dir, region_to_test, Sys.Date()-1), row.names=F)

  R0_val = Rtoutputs$R0
  plotout_noage = bind_rows(plotout_noage, Rt)
  all_times = as.numeric(unique(plotout_noage$Date) - as.Date('2020-01-14'))

  icuf = get_time_varying_pars(new_pars[grep(names(new_pars), pattern="ICU_values_*")],
      new_pars[['ICU_max']],
      new_pars[['ICU_min']],
      new_pars[grep(names(new_pars), pattern="ICU_changepoint_values_*")],
      all_times)

  icu_frac_df = data.frame(Date = as.Date('2020-01-14') + all_times, Compartment='ICU fraction', Cases=icuf, Source='Simulation')

  plotout_noage = bind_rows(plotout_noage, icu_frac_df)

  init_infect = round(popfinal[[region_to_test]] * new_pars[[paste0('num_init_', 1)]])

  print("Plotting")
  ggplot(plotout_noage, aes(x=Date, y=Cases, color=Source, fill=Source)) +
     stat_summary(fun=function(z){quantile(z,0.5,type=3)}, geom="line", size=0.75) +
     stat_summary(fun.min=function(z){quantile(z,0.025, type=3)}, fun.max=function(z){quantile(z,0.975, type=3)}, geom="ribbon", alpha=0.3, color=NA) +
     geom_point(plotting_data, mapping=aes(x=Date, y=Cases, fill=Source), size=0.9, alpha=0.7, pch=21) +
     facet_wrap(~factor(Compartment, levels=plot_order), scales='free_y', ncol=3) +
     theme_bw() +
     ylab('') +
     scale_x_date(date_breaks = "months", date_labels="%b") +
     ggtitle(sprintf("Region %s\nR0: %s\nInitial prevalence: %s\nloglik: %s", 
      region_to_test, R0_val, init_infect, loglik)) +
     scale_color_brewer(palette='Set1', limits=compartment_order) +
     scale_fill_brewer(palette='Set1', limits=compartment_order)
  ggsave(sprintf('%s/projection_%s_%s.png', human_output_dir, model_name, region_to_test), width=14, height=10)

  ggplot(plotout_noage %>% filter(Date >= project_zoom), aes(x=Date, y=Cases, color=Source, fill=Source)) +
     stat_summary(fun=function(z){quantile(z,0.5,type=3)}, geom="line", size=0.75) +
     stat_summary(fun.min=function(z){quantile(z,0.025, type=3)}, fun.max=function(z){quantile(z,0.975, type=3)}, geom="ribbon", alpha=0.3, color=NA) +
     geom_point(plotting_data%>% filter(Date >= project_zoom), mapping=aes(x=Date, y=Cases, fill=Source), size=0.9, alpha=0.7, pch=21) +
     facet_wrap(~factor(Compartment, levels=plot_order), scales='free_y', ncol=3) +
     theme_bw() +
     ylab('') +
     scale_x_date(date_breaks = "months", date_labels="%b") +
     ggtitle(sprintf("Region %s\nR0: %s\nInitial prevalence: %s\nloglik: %s", 
      region_to_test, R0_val, init_infect, loglik)) +
     scale_color_brewer(palette='Set1', limits=compartment_order) +
     scale_fill_brewer(palette='Set1', limits=compartment_order)
  ggsave(sprintf('%s/projection_zoom_%s_%s.png', human_output_dir, model_name, region_to_test), width=14, height=10)
  }

get_time_varying_pars = function(cpoint_logit_values, logit_max,logit_min, cpoint_times, all_times){
    cpoint_values = mapply(FUN=unlogit, cpoint_logit_values, pmax=logit_max, pmin=logit_min)

    df = data.frame(time=cpoint_times, values=cpoint_values)
    df = full_join(df, data.frame(time=all_times), by='time') %>% 
        arrange(time)
    na.approx(df) %>%
        as.data.frame() %>%
        mutate(values = case_when(time < min(cpoint_times) ~ cpoint_values[1],
                                  time > max(cpoint_times) ~ cpoint_values[length(cpoint_values)],
                                  T ~ values)) %>%
        arrange(time) %>%
        select(values) %>%
        unlist()
}
