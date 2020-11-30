
run_epinow2 <- function(dat_df,  # Data used in estimation
                        obs_colname, # Name of column holding observations within data_df
                        dat_type, # Can be 'cases', 'deaths' or 'hospitalizations'
                        prior_smoothing_window = 7, # Smoothing window to use on the prior. Default is 7
                        dbug, # If true, run really short chains.
                        #midway = FALSE,
                        output_folder = 'rough-rt-approach'){
  
  
  if(!dir.exists(output_folder)){dir.create(output_folder)}
  if(!dir.exists(paste0(output_folder, '/figs'))){dir.create(paste0(output_folder, '/figs'))}
  
  ## Set delay distributions for input into EpiNow2 -------------------------------------------
  incubation_period <- EpiNow2::get_incubation_period('SARS-CoV-2', 'lauer')
  incubation_period$notes = 'lauer et al'
  # See EpiNow2::incubation_periods for source info
  
  generation_time <- EpiNow2::get_generation_time('SARS-CoV-2', source = 'ganyani')
  generation_time$notes = 'ganyani et al'
  ## See EpiNow2::generation_intervals for source info
  
  cd_post <- read_rds('../data/fitted_delays/delay_symptom_to_report_posterior_lognorm.rds') %>%
    bind_rows(.id = 'par')
  case_rep_delay <- list(mean = median(cd_post$mu), ## Very rough estimates based on carline 
                         mean_sd = sd(cd_post$mu), 
                         sd = median(cd_post$sigma),
                         sd_sd = sd(cd_post$sigma),
                         max = 30,
                         notes = 'Phils fitted estimates')
  
  death_rep_delay <- list(mean = 2.86, ## From Linton et al. Table 2 (J. Clin. Med. 2020, 9, 538; doi:10.3390/jcm9020538)
                          mean_sd = 1.98,   
                          sd = 0.53,
                          sd_sd = 2.03,
                          max = 30,
                          notes = 'From Linton et al. Table 2. (J. Clin. Med. 2020, 9, 538; doi:10.3390/jcm9020538)')
  
  hospital_rep_delay <- list(mean = 2.21, ## Fitted to HK data, similar to Linton et al
                             mean_sd = 0.1, 
                             sd = 0.49,
                             sd_sd = 0.1,
                             max = 30,
                             notes = 'From HK data. Results were similar to Linton et al.')
  
  stopifnot(dat_type %in% c('cases', 'deaths', 'hospitalizations'))
  if(dat_type == 'cases')   delay <- case_rep_delay  
  if(dat_type == 'deaths')   delay <- death_rep_delay 
  if(dat_type == 'hospitalizations')   delay <- hospital_rep_delay 
  
  
  write_rds(generation_time, sprintf('%s/gen_interval.rds', output_folder))
  write_rds(incubation_period, sprintf('%s/incubation_pd.rds', output_folder))
  write_rds(delay, sprintf('%s/delay.rds', output_folder))
  
  
  # Plot the specified distributions
  png(sprintf('%s/figs/specified_distributions.png', output_folder), width = 7, height = 7, units = 'in', res = 300)
  par(mfrow = c(2,2))
  xx = seq(0, 30, by = 0.01)
  ## Reporting delay --
  plot(xx, dlnorm(xx, delay$mean, delay$sd), 
       type = 'l', main = sprintf('Rep. delay (lognormal) logmean=%2.2f, logsd=%2.2f\nsource - %s', delay$mean, delay$sd, delay$notes), 
       xlab = 'days', ylab = 'dens')
  ## Generation interval ---
  plot(xx, 
       dgamma(xx, shape = with(generation_time, get_shape(mean, sd^2)), 
              rate = with(generation_time, get_rate(mean, sd^2))), 
       type = 'l', main = sprintf('Gen int. (gamma) mean=%2.2f, sd=%2.2f\nsource - %s', generation_time$mean, generation_time$sd, generation_time$notes),  xlab = 'days', ylab = 'dens')
  ## Incubation period --
  plot(xx, 
       dlnorm(xx, incubation_period$mean, incubation_period$sd), 
       type = 'l', main = sprintf('Incubation (lognormal) mean=%2.2f, sd=%2.2f\nsource - %s', incubation_period$mean, incubation_period$sd, incubation_period$notes), xlab = 'days', ylab = 'dens')
  dev.off()
  
  
  ## Input into inference model  -------------------------------------------
  ## Write a wrapper to reformat the desired synthetic data for input into epiEstim
  format_dat <- function(obs_colname, odf = dat_df){
    odf[,c('date', obs_colname)] %>%
      setNames(c('date', 'confirm')) %>%
      data.table::as.data.table() %>%
      return()
  }
  
  
  ## Fit to synthetic case observations
  rt_estimates <- EpiNow2::epinow(reported_cases = format_dat(obs_colname, dat_df), 
                                  generation_time = generation_time,
                                  delays = list(reporting_delay = delay,
                                                incubation_period = incubation_period), 
                                  method = 'exact',
                                  CrIs = c(.8, .9, .95),
                                  prior_smoothing_window = prior_smoothing_window,
                                  rt_prior = list(mean = 2, sd = 1), 
                                  horizon = 0,
                                  samples = ifelse(dbug, 10, 2000), 
                                  stan_args = list(warmup = ifelse(dbug, 70, 500), 
                                                   control = list(adapt_delta = 0.9),
                                                   cores = 4),
                                  target_folder = paste0(output_folder))
  
  
  # if(!debug & midway){
  #   ## Fit to synthetic case observations
  #   rt_estimates <- EpiNow2::epinow(reported_cases = format_dat(obs_colname, dat_df), 
  #                                   generation_time = generation_time,
  #                                   reporting_delay = delay,
  #                                   incubation_period = incubation_period, 
  #                                  # prior_smoothing_window = prior_smoothing_window,
  #                                   rt_prior = list(mean = 2, sd = 1), horizon = 0,
  #                                   samples = 2000, warmup = 500, cores = 4,
  #                                   chains = 4, verbose = TRUE,
  #                                   target_folder = paste0(output_folder))
  # }
}
