fit_lnorm_quantiles <- function(
  quantile.values= c(0, 2, 5,8,20),
  quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975),
  mean.delay=5.6){
  # Fit a lognormal to quantile data on time from symptom onset to hospital admission
  # from Belgium (Faes, et al.)
  fun_to_optimize = function(par){
      sd.est = par
      mu.est = log(mean.delay) - sd.est^2/2
      q = qlnorm(quantiles, meanlog=mu.est, sdlog = sd.est)
      return (-sum((q - quantile.values)^2))
  }
  best.sd = optimize(f=fun_to_optimize, 
           lower=0.1,
           upper=3,
           maximum = T)$maximum
  best.mean = log(mean.delay) - best.sd^2/2

  return(list(meanlog.fit=best.mean, sdlog.fit=best.sd))

}

run_epinow2 <- function(dat_df,  # Data used in estimation
                        obs_colname, # Name of column holding observations within data_df
                        prior_smoothing_window = 7, # Smoothing window to use on the prior. Default is 7
                        dbug, # If true, run really short chains.
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

  best_fit = fit_lnorm_quantiles()

  # DELAY FROM SYMPTOM ONSET TO HOSPITAL ADMISSION
  delay <- list(mean =best_fit[['meanlog.fit']],
                         mean_sd = 0.1, 
                         sd = best_fit[['sdlog.fit']],
                         sd_sd = 0.1,
                         max = 15,
                         notes = 'Fitted to Faes et al.')

  write_rds(generation_time, sprintf('%s/gen_interval.rds', output_folder))
  write_rds(incubation_period, sprintf('%s/incubation_pd.rds', output_folder))
  write_rds(delay, sprintf('%s/delay.rds', output_folder))
  
  
  # Plot the specified distributions
  png(sprintf('%s/figs/specified_distributions.png', output_folder), width = 7, height = 7, units = 'in', res = 300)
  par(mfrow = c(2,2))
  xx = seq(0, 30, by = 0.01)
  ## Reporting delay --
  if (length(delay) > 0){
      plot(xx, dlnorm(xx, delay$mean, delay$sd), 
       type = 'l', main = sprintf('Rep. delay (lognormal) logmean=%2.2f, logsd=%2.2f\nsource - %s', delay$mean, delay$sd, delay$notes), 
       xlab = 'days', ylab = 'dens')
  }

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
                                  delays = list(reporting_delay= delay, incubation_period=incubation_period), 
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
}
