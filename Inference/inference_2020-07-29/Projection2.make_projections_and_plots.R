library(dplyr)
library(tidyr)
library(pomp)
library(rmutil)
library(foreach)
library(gridExtra)
library(digest)
library(ggplot2)
library(stringr)
library(data.table)


select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

root <- '../../'
source(file.path(root, '_covid_root.R'))

covid_set_root(root)
default_par_file = './final_mle_pars.csv'
fit_df = read.csv('./final_points.csv')


source("./input_file_specification.R") # script to read in necessary files
output_dir = projection_dir
source(covid_get_path(inference_file))
source("./set_up_covariates_and_data.R")
source(covid_get_path(projection_functions))


# Projection specifications: all projections
reference_date <- as.Date(t_ref) # All numeric times in reference to this date
end_projection_date <- Sys.Date()
start_projection_date <- simstart + reference_date

R0_values = read.csv(paste0(projection_dir, 'R0_values.csv')) %>%
  mutate(Date=reference_date+time, Region=as.character(Region))

deltaT = 0.1 # timestep for projections  

regional_aggregation = F # If true, do statewide estimates
initialize=T
region_order = c('northcentral','central','northeast','southern', 'chicago')
civis_data = read.csv(data_filename)  %>% 
    mutate(Date=as.Date(date), total_deaths = hosp_deaths+nonhosp_deaths) %>%
    filter(Date >= as.Date('2020-03-16'))

get_idph = function(date, region){
    idph[which(idph$Date==date & idph$restore_region==region), 'total_deaths']
}

idph = read.csv(idph_filename) %>% 
    mutate(Date=as.Date(test_date), restore_region=region, total_deaths=inc_deaths) %>%
    select(Date, restore_region, total_deaths)

civis_data = civis_data %>% mutate(
    source=case_when((is.na(total_deaths)) ~ 'Public',
                                (!is.na(total_deaths))~ 'IDPH line list'),
    total_deaths=case_when((is.na(total_deaths)) ~ as.numeric(mapply(get_idph, Date, restore_region)),
                                (!is.na(total_deaths))~ as.numeric(total_deaths)),
    incident_hospitalizations=new_ll_hospitalizations,
    confirmed_covid_nonicu=covid_non_icu
                                )
IL_civis = civis_data %>%
    filter(restore_region != 'chicago') %>%
    group_by(date, Date, source, source_hosp_deaths) %>%
    summarize(emr_deaths = sum(emr_deaths),
              nonhosp_deaths = sum(nonhosp_deaths),
              confirmed_covid_icu = sum(confirmed_covid_icu),
              total_deaths=sum(total_deaths),
              incident_hospitalizations=sum(new_ll_hospitalizations),
              confirmed_covid_nonicu=sum(covid_non_icu)) %>%
        ungroup() %>%
    mutate(restore_region = 'Illinois')
    
civis_data = bind_rows(civis_data, IL_civis) %>%
    mutate(confirmed_covid_nonicu = if_else(Date < as.Date('2020-05-06'), NA_integer_, confirmed_covid_nonicu))

pop_totals = as.numeric(lapply(FUN=sum, population_list))
chicago_pop = pop_totals[5]
pop_totals = pop_totals[1:4]
pop_totals = c(pop_totals, chicago_pop, sum(pop_totals))
names(pop_totals) = c(region_order[1:4], 'chicago', 'Illinois')

## Simulate up to current date

# Loop through and do simulations for each parameter set
for (parset_id in 1:max(fit_df$parset)){
 
    pars_to_use = fit_df %>% filter(parset == parset_id) %>% select(-X) %>% unlist(use.names=T)
    pars[names(pars_to_use)] = pars_to_use
    num_sims = pars$num_sims

    end_projection_date <- as.Date(final_projection_date) #Sys.Date()
    start_projection_date <- simstart + reference_date

    ## Set up simulation timng params
    t0_sim = as.numeric(start_projection_date - reference_date)
    simend = as.numeric(end_projection_date - reference_date)
    pars$tmin=t0_sim
    pars$tmax = simend 

    sim_full = simulate_pomp_covid(
      n_regions = pars$n_regions,
      n_age_groups = n_age_groups,
      nsim=num_sims, 
      input_params = pars,
      delta_t = deltaT,
      population=population_list,
      beta_scales=beta_scales,
      beta_covar=beta_covariate_column,
      contacts=pomp_contacts,
      frac_underreported=fraction_underreported,
      rprocess_Csnippet = rprocess_snippet,
      rinit_Csnippet = rinit_snippet,
      rmeasure_Csnippet=rmeasure_snippet,
      initialize=T,
      obsnames =observed_names
    )

######################################## Bug in reinitialization...
#    print("Reinitializing")
    # Get new ICs for new starting point
#    param_values = redistribute_hospitalizations()
    # Simulate from different start and end projection date with new ICs
#    start_projection_date <- Sys.Date()
#    end_projection_date <- as.Date(final_projection_date)
    
    ## Set up simulation timng params
#    t0_sim = as.numeric(start_projection_date - reference_date)
#    simend = as.numeric(end_projection_date - reference_date)
#    pars$tmin=t0_sim
#    pars$tmax = simend 

#    init_params = param_values

    
#    print("Simulating from new ICs")
#    # Simulate once for every simulation in first half
#    foreach(sim=1:nrow(param_values), .combine='rbind') %do% {
#        init_params = as.numeric(param_values[sim, ])
#        names(init_params) = paste0(names(param_values), '.0')
        
#        temp_sim = simulate_pomp_covid(
#            n_age_groups=nrow(population_list[[1]]),
#            n_regions = pars$n_regions,
#            nsim=1, 
#            input_params=pars,
#            delta_t=deltaT,
#            population=population_list,
#            contacts=pomp_contacts,
#            beta_scales=beta_scales,
#            beta_covar=beta_covariate_column,
#            frac_underreported = fraction_underreported,
#            rprocess_Csnippet = rprocess_snippet,
#            rinit_Csnippet = rinit_snippet,
#            rmeasure_Csnippet=rmeasure_snippet,
#            obsnames=observed_names,
#            initialize = F,
#            initial_params=init_params

#        )
#        temp_sim$raw_simulation_output$.id = sim
#        temp_sim$raw_simulation_output
#    } -> sim_raw2
    
    # Fix accumulators
#        param_values = sim_full$raw_simulation_output %>% filter(time == max(time)) %>% 
#            select(-time, -.id)
        # get initial_value for all accumulators
#        accum = c(names(param_values)[grep('^new_', names(param_values))], 
#                  names(param_values)[grep('^Inc_', names(param_values))],
#                  names(param_values)[grep('^Obs', names(param_values))])
#        accum_init = param_values[accum]
#        t0 = min(sim_raw2$time)
#        r_true = 1
#        for (r in 1:nrow(sim_raw2)){
#            if(sim_raw2[r, 'time'] == t0){
#               sim_raw2[r, accum] = accum_init[r_true, accum]
#               r_true = r_true + 1
#            }
        
#        }
    
#    print("Binding results")
    # Bind results together, don't forget to add in a parset ID
#        sim_full$raw_simulation_output$.id = as.numeric(sim_full$raw_simulation_output$.id)
#        sim_full$raw_simulation_output = bind_rows(sim_full$raw_simulation_output, sim_raw2)
####################################################################

        sim_full$raw_simulation_output$.id = as.numeric(sim_full$raw_simulation_output$.id)
        alloutputs = process_pomp_covid_output(sim_full, agg_regions=F)
        
        plotout = alloutputs$plotting_output %>% ungroup() %>%    
            mutate(Date=reference_date + Time)
        df_infections_summary <- plotout %>%
            group_by(Time, Compartment, Region) %>%
            filter(Compartment %in% c('E', 'A', 'P', 'IS', 'IM', 'IH', 'IC')) %>%
            group_by(SimID, Date, Region) %>%
            arrange(Time) %>%
            summarise(Compartment = 'Prevalence',
                      Cases = sum(Cases)) %>%
          ungroup()
        
        chicago_out = bind_rows(plotout, df_infections_summary) %>% 
            filter(Region==5) %>%
            mutate(restore_region=region_order[as.numeric(Region)]) %>%
            select(-Time)
        
        plotout = bind_rows(plotout, df_infections_summary) %>% 
            filter(Region!=5) %>% # Remove Chicago
            mutate(restore_region=region_order[as.numeric(Region)]) %>%
            select(-Time)
        
        statewide = plotout %>%
            group_by(Date, SimID, Compartment) %>%
            summarize(Cases=sum(Cases)) %>%
            mutate(Region='0', restore_region='Illinois') %>%
            ungroup()
        
        plotout = bind_rows(plotout, statewide, chicago_out) %>% filter(!is.na(Compartment))
        plotout$parset = parset_id
        
        Rtdf = add_R0_to_output(plotout, R0_values, parset)
        plotout = rbind(plotout, Rtdf)
        
        write.csv(plotout, paste0(model_name, '_parset_', parset_id, '.csv'), row.names=F)
        rm(sim_full)
        
    }

file_list <- list.files(pattern = paste0(model_name, '_parset_.*.csv'))
cons <- rbindlist(lapply(file_list, fread))
cons <- as.data.frame(cons)
write.csv(cons, paste0(output_dir,model_name, '_all.csv'))    


### Plotting projections ###
outfile= paste0(output_dir, model_name, '_all.fitplot.png')
outfile2 = paste0(output_dir, model_name, '_all.states.png')

region_plot_order=c('northeast','southern','central','northcentral', 'chicago', 'Illinois')

foreach(rnum=1:6) %do%{
    region_to_plot = region_plot_order[rnum]

    civisplot = civis_data %>% 
        mutate(deaths=total_deaths) %>%
        gather(Compartment, Cases, deaths, confirmed_covid_icu, confirmed_covid_nonicu, emr_deaths, nonhosp_deaths, incident_hospitalizations) %>% 
        filter(restore_region==region_to_plot) %>%
        mutate(source=case_when((Compartment=='confirmed_covid_icu') ~ 'emresource',
                                (Compartment=='deaths') ~ source,
                                (Compartment=='emr_deaths') ~ 'emresource',
                                (Compartment=='nonhosp_deaths') ~ 'IDPH line list',
                                (Compartment=='incident_hospitalizations') ~ 'IDPH line list',
                                (Compartment=='confirmed_covid_nonicu') ~ 'emresource')) %>%
        mutate(Compartment=case_when((Compartment=='deaths')~'All reported deaths',
                                     (Compartment=='confirmed_covid_icu')~'Confirmed ICU cases',
                                     (Compartment=='nonhosp_deaths')~'Non-hospital deaths',
                                     (Compartment=='emr_deaths')~'Hospitalized deaths',
                                     (Compartment=='incident_hospitalizations') ~ 'Incident hospitalizations',
                                     (Compartment=='confirmed_covid_nonicu') ~ 'Confirmed non-ICU hospital cases' ))

    plotout %>% 
        filter(Compartment %in% c('Reported deaths', 'ObsICU', 'nNHD', 'ObsHospDeaths', 'new_hospitalizations', 'ObsHosp')) %>% 
        filter(Date<as.Date('2020-09-30'), restore_region==region_to_plot) %>%
        mutate(Compartment=case_when((Compartment=='Reported deaths')~'All reported deaths',
                                     (Compartment=='ObsICU')~'Confirmed ICU cases',
                                     (Compartment=='nNHD')~'Non-hospital deaths',
                                     (Compartment=='ObsHospDeaths')~'Hospitalized deaths',
                                     (Compartment=='new_hospitalizations') ~ 'Incident hospitalizations',
                                     (Compartment=='ObsHosp') ~ 'Confirmed non-ICU hospital cases')) %>%
        ggplot(aes(x=Date, y=Cases)) + 
        stat_summary(fun.ymin=function(z){quantile(z,0.025)}, fun.ymax=function(z){quantile(z,0.975)}, geom="ribbon", color="black", alpha=0.1) +
        stat_summary(fun.y=median, geom="line", color="black", size=1) +
        geom_point(civisplot, mapping=aes(x=Date, y=Cases, fill=source), color='black', pch=21, size=3, alpha=0.5) +
        geom_smooth(civisplot, mapping=aes(x=Date, y=Cases), color='firebrick', se=F, span=0.75) +
        facet_wrap(~Compartment, scales='free_y', ncol=6) +
        ylab(region_to_plot) + theme_bw() +
        scale_fill_brewer(palette='Dark2') -> p
    p
} -> plots

g = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=6)
ggsave(outfile, g, width=16, height=16)


region_plot_order=c('northeast','southern','central','northcentral', 'chicago', 'Illinois')
order_to_plot_states = c('Prevalence',
                         'Incidence',
                         'Recovered',
                         'Susceptible',
                         'Rt',
                         'Incident deaths',
                         'Non-ICU hospital beds occupied',
                         'ICU beds occupied')

foreach(rnum=1:6) %do%{
    region_to_plot = region_plot_order[rnum]
    popregion = pop_totals[region_to_plot]
plotout %>% 
    filter(Compartment %in% c('Incidence', 'Prevalence','IH','IC', 'nD', 'R', 'S', 'Rt')) %>% 
    mutate(Compartment = case_when((Compartment=='Incidence') ~ 'Incidence',
                                   (Compartment=='Prevalence') ~ 'Prevalence',
                                   (Compartment == 'IH') ~ 'Non-ICU hospital beds occupied',
                                   (Compartment == 'IC') ~ 'ICU beds occupied',
                                   (Compartment == 'nD') ~ 'Incident deaths',
                                   (Compartment == 'R') ~ 'Recovered',
                                   (Compartment == 'S') ~ 'Susceptible',
                                   (Compartment == 'Rt') ~ 'Rt')) %>%
    filter(Date<as.Date('2020-11-30'), restore_region==region_to_plot) %>%
    mutate(Cases = case_when((Compartment %in% c('Prevalence','Recovered', 'Susceptible')) ~ Cases / popregion,
                            (!Compartment %in% c('Prevalence','Recovered', 'Susceptible')) ~ Cases)) %>%
    mutate(Compartment = factor(Compartment, levels=order_to_plot_states))  -> incidence_plots

incidence_plots %>%
    filter(Compartment!='Susceptible') %>%
    ggplot(aes(x=Date, y=Cases)) + 
    stat_summary(fun.ymin=function(z){quantile(z,0.025)}, fun.ymax=function(z){quantile(z,0.975)}, geom="ribbon", color="black", alpha=0.1) +
    stat_summary(fun.y=median, geom="line", color="black", size=1) +
    facet_wrap(~Compartment, scales='free_y', ncol=8) +
    ylab(region_to_plot) + theme_bw()->p
p
} -> plots

g = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=6)
ggsave(outfile2, g, width=32, height=20)


### Format for civis ###

regions = c('Illinois','restore_northcentral', 'restore_central', 'restore_northeast', 'restore_southern', 'chicago')

covid_set_root(root)

intervention = model_name

agg = plotout %>% select(Date, Region, Compartment, Cases)
agg$Region = sapply(agg$Region, FUN=function(x){regions[as.numeric(x) + 1]})
agg$intervention=intervention

ptilehigh <- function(x){return(quantile(x, 0.975))}
ptilelow <- function(x){return(quantile(x, 0.025))}
ptilemed <- function(x){return(quantile(x, 0.5))}


Rt_all = agg %>% filter(Compartment=='Rt') %>%
    group_by(Date, intervention, Region) %>%
    summarize(Rt_median=ptilemed(Cases),
          Rt_lower=ptilelow(Cases),
          Rt_upper=ptilehigh(Cases))

prev = agg %>% filter(Compartment=='Prevalence') %>%
    group_by(Date, intervention, Region) %>%
    summarize(cases_median=ptilemed(Cases),
          cases_lower=ptilelow(Cases),
          cases_upper=ptilehigh(Cases))

inc = agg %>% filter(Compartment=='Incidence') %>%
    group_by(Date, intervention, Region) %>%
    summarize(cases_new_median=ptilemed(Cases),
          cases_new_lower=ptilelow(Cases),
          cases_new_upper=ptilehigh(Cases))

deaths = agg %>% filter(Compartment=='nD') %>%
    group_by(Date, intervention, Region) %>%
    summarize(deaths_median=ptilemed(Cases),
           deaths_lower=ptilelow(Cases),
          deaths_upper=ptilehigh(Cases))
    
deaths_det = agg %>% filter(Compartment=='Reported deaths') %>%
    #group_by(SimID, Date, intervention, Region) %>%
    #summarize(Cases = sum(Cases)) %>%
    mutate(Compartment = 'deaths_det') %>%
    #ungroup() %>%
    group_by(Date, intervention, Region) %>%
    summarize(deaths_det_median=ptilemed(Cases),
           deaths_det_lower=ptilelow(Cases),
          deaths_det_upper=ptilehigh(Cases))


hosp = agg %>% filter(Compartment=='IH') %>%
    group_by(Date, intervention, Region) %>%
    summarize(hosp_bed_median=ptilemed(Cases),
          hosp_bed_lower=ptilelow(Cases),
          hosp_bed_upper=ptilehigh(Cases))

icu = agg %>% filter(Compartment=='IC') %>%
    group_by(Date, intervention, Region) %>%
    summarize(icu_median=ptilemed(Cases),
          icu_lower=ptilelow(Cases),
          icu_upper=ptilehigh(Cases))

vent = icu %>% 
    mutate(vent_median=0.74*icu_median,
           vent_lower=0.74*icu_lower,
           vent_upper=0.74*icu_upper) %>%
    select(-icu_median, -icu_lower, -icu_upper)

recov = agg %>% filter(Compartment=='R') %>%
    group_by(Date, intervention, Region) %>%
    summarize(recovered_median=ptilemed(Cases),
          recovered_lower=ptilelow(Cases),
          recovered_upper=ptilehigh(Cases))

final_output = inner_join(prev, inc, by=c('Date', 'intervention', 'Region')) %>% 
    inner_join(deaths, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(deaths_det, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(hosp, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(icu, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(vent,by=c('Date', 'intervention', 'Region')) %>%
    inner_join(recov, by=c('Date', 'intervention', 'Region')) %>%
    inner_join(Rt_all, by=c('Date', 'intervention', 'Region')) %>%
    filter(Date >= as.Date('2020-03-13')) %>% 
    arrange(Region, intervention, Date)

names(final_output) = c('date', 'scenario_name', 'geography_modeled', names(final_output[4:ncol(final_output)]))
write.csv(final_output, sprintf('%s/UChicago_%s_%s.csv',output_dir, intervention, Sys.Date()), row.names=F)

plotout %>% 
    filter(Compartment %in% c('D','R'), Date==Sys.Date()) %>% 
    group_by(restore_region, Compartment, Date) %>%
    summarize(mean_count=mean(Cases)) ->ifr_df

write.csv(ifr_df, sprintf('%s/ifr_%s.csv',output_dir, Sys.Date()), row.names=F)
