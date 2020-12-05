covid_source('utils.R')

library(stringr)
#' Simulate COVID pomp model
simulate_pomp_covid <- function(
  n_regions,
  input_params, 
  delta_t, 
  emr_report_covar,
  icu_covars,
  inc_hosp,
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
  params <- simulate_pomp_covid__init_parameters(input_params)

  # Set up covariate table for pomp
  covar_table_interventions <- simulate_pomp_covid__init_covariate_table(
  input_params,
  icu_covars,
  emr_report_covar,
  inc_hosp
  )

  covar_table <- covariate_table(
    covar_table_interventions,
    order = "constant",
    times = "time"
  )

  # Actually call pomp
  state_names <- simulate_pomp_covid__init_state_names(n_age_groups, n_regions, subcompartment_df)
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
      statenames = c(state_names),
      paramnames = names(params),
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
      state = c("preIC_", "IC_"),
      index = c(alpha_preIC,alpha_IC

      )
    )
  })
}

simulate_pomp_covid__init_state_names <- function(n_age_groups, n_regions, subcompartment_df) {
  state_names = c()
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
  input_params
) {
  with(input_params, {
    # Set up parameters for model
    params = c(
      # Control params
      "n_regions" = n_regions,
      "region_to_test" = region_to_test,
      "alpha_IC" = alpha_IC, 
      "alpha_preIC" = alpha_preIC
    )
    params
  })
}

simulate_pomp_covid__init_covariate_table <- function(input_params, 
  icu_covars,
  emr_report_covar,
  inc_hosp) {
    
  covar_table_interventions <- data.frame(time = c(1:input_params$tmax))
  
  # add emr_reporting
  covar_table_interventions = left_join(covar_table_interventions, emr_report_covar, by='time')
  covar_table_interventions = left_join(covar_table_interventions, icu_covars, by='time')
  covar_table_interventions = left_join(covar_table_interventions, inc_hosp, by='time')
  covar_table_interventions %>% mutate(int_time=time)
}
 
process_pomp_covid_output <- function(sim_result, agg_regions=T) {
  num_age_classes = 1
  num_regions = 11
  max_subcompartments = 3
  subcompartments = c("IC", "preIC")
  obs_prefix = c("ObsDeaths", "ObsHosp", "ObsHospDeaths")

  sub_suffixes = expand.grid(subcompartments, 1:max_subcompartments, 1:num_age_classes, 1:num_regions)
  sub_suffixes$state_name = paste(sub_suffixes$Var1, sub_suffixes$Var2, sub_suffixes$Var3, sub_suffixes$Var4,  sep='_')
  sub_suffixes$Var1 = gsub('[[:digit:]]+', '', sub_suffixes$Var1)

  obs_suffixes = expand.grid(obs_prefix, 1:num_regions)
  obs_suffixes$state_name = paste(obs_suffixes$Var1, obs_suffixes$Var2,  sep='_')

  name_map = c(as.vector(sub_suffixes$Var1), as.vector(obs_suffixes$Var1))
  names(name_map) = c(sub_suffixes$state_name, obs_suffixes$state_name)

  region_map = c(as.vector(sub_suffixes$Var4), as.vector(obs_suffixes$Var2))
  names(region_map) = names(name_map)

  age_map = c(as.vector(sub_suffixes$Var3), as.vector(rep(-1, length(as.vector(obs_suffixes$Var2)))))
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