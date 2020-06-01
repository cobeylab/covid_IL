covid_source('utils.R')
covid_source('Forecasting/simulation_statewide.R')

#' Simulate COVID pomp model
make_pomp_object_covid <- function(
  n_age_groups,
  n_regions,
  input_params, 
  delta_t, 
  contacts, 
  population_list, 
  nu_scales,
  beta_scales,
  frac_underreported,
  fitstart,
  dmeasure_Csnippet,
  rprocess_Csnippet,
  data,
  time_column,
  transformations,
  rinit_Csnippet,
  inherit_parameters=TRUE,
  format = 'data.frame',
  obsnames = NULL,
  inherited_params = NULL
) {
  library(pomp)
  print('Making subcompartments')
  # Initialize states
  subcompartment_df <- simulate_pomp_covid__init_subcompartment_df(input_params)
  # Intervention df
  intervention_df <- simulate_pomp_covid__init_intervention_df(input_params)
  n_interventions <- nrow(intervention_df)
  
  if (inherit_parameters)
  {
    params <- inherited_params

  } else{
    # Set up parameters for pomp
    params <- simulate_pomp_covid__init_parameters(
      n_age_groups, n_interventions, input_params, contacts, population_list
    )
    print('params initialized')
  }

  # Set up covariate table for pomp
  covar_table_interventions <- simulate_pomp_covid__init_covariate_table(
    input_params, intervention_df, nu_scales, beta_scales, frac_underreported
  )

  covar_table <- covariate_table(
    covar_table_interventions,
    order = "constant",
    times = "time"
  )

print('covars initialized')
  # Actually call pomp
  model_path <- covid_get_path('POMP/model_scripts/statewide')
 
  ## Set up state variables by age group and then by regions 
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


print('calling pomp')
  pomp_object <- pomp(
    times = time_column,
    t0 = fitstart,
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
    obsnames = obsnames)
  pomp_object
}

get_civis_data <- function(civis_filename){
  data = read.csv(civis_filename)
  data = data %>% select(date, new_deaths, confirmed_covid_icu) %>% filter(!is.na(new_deaths))
  names(data) = c('date','ObsDeaths', 'ObsICU')
  data$ObsDeaths = round(data$ObsDeaths)
  data$ObsCases = NA
  data$time = as.numeric(as.Date(data$date) - as.Date('2020-01-14'))
  data = rbind(c(NA, NA, NA, NA, min(data$time)-1), data)
  data
}


run_pfilter_and_output_result <- function(n_reps_pfilter, 
  n_particles_pfilter, 
  jobid,
  arrayid,
  output_dir,
  pomp_object,
  dmeasFile,
  simstart,
  min_data_time,
  intervention_start,
  calculate_mean=T){
  require(foreach)

  chainID = sprintf("%s_%s", jobid, arrayid)
  output_filename = sprintf("%s/output_pfilter_%s_%s.csv", output_dir, jobid, arrayid)
  t = Sys.time()

  if (calculate_mean){
    ll <- logmeanexp(replicate(n_reps_pfilter,logLik(pfilter(pomp_object,Np=n_particles_pfilter))),se=TRUE)
    runtime = Sys.time() - t

    finaloutput = (data.frame(as.list(pomp_object@params),
                        loglik=ll[1],
                        loglik_se=ll[2],
                        chain = chainID,
                        model=dmeasFile,
                        t0=simstart,
                        data_start=min_data_time,
                        intervention_start=intervention_start,
                        runtime=runtime,
                        nparticles=n_particles_pfilter,
                        n_reps=n_reps_pfilter))    
    } else{
    foreach (i=1:n_reps_pfilter, .combine='rbind') %do%{
      ll <- logLik(pfilter(pomp_object,Np=n_particles_pfilter))
      runtime = Sys.time() - t
      data.frame(as.list(pomp_object@params),
        loglik=ll,
        loglik_se=NA,
        chain = chainID,
        model=dmeasFile,
        t0=simstart,
        data_start=min_data_time,
        intervention_start=intervention_start,
        runtime=runtime,
        nparticles=n_particles_pfilter,
        n_reps=1)
      } -> finaloutput
    }
  write.table(finaloutput, file = output_filename, sep = ",",col.names = TRUE, row.names=FALSE, append=TRUE)
}

run_mif_and_output_result <- function(
  jobid,
  arrayid,
  output_dir,
  pomp_object,
  dmeasFile,
  simstart,
  min_data_time,
  intervention_start,
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
  output_filename = sprintf("%s/output_mif_%s_%s.csv", output_dir, jobid, arrayid)

  saveRDS(mf, paste0(output_filename, '.rds'))
}

convert_date_to_time = function(d1, dref){
  as.numeric(as.Date(d1) - as.Date(dref))
}