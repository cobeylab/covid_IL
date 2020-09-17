redistribute_hospitalizations = function(){
    # Prepare initial states
    param_values = sim_full$raw_simulation_output %>% filter(time == max(time)) %>% 
        select(-time, -.id) %>% 
        select(-observed_names) # remove observed parameters
    
    final_param_names = paste0(names(param_values), '.0')
 
    # Re-initialize ICs

    for (sim in 1:nrow(param_values)){
        # For each simulation, pick an IFR
        ifr = runif(1, 0.0065, 0.009)
        for (region in 1:5){
            regex= paste0('^R_[1-9]_', region)
            param_names = names(param_values)[grep(regex, names(param_values))]
            
            regex_S = paste0('^S_[1-9]_', region)
            S_names = names(param_values)[grep(regex_S, names(param_values))]
            
            regex_D = paste0('^D_[1-9]_', region)
            D_names = names(param_values)[grep(regex_D, names(param_values))]
            
            latent_deaths = sum(param_values[sim, D_names])
            region_name = region_order[region]

            # Something is broken about this part that I can't figure out
            #new_recovered = round(latent_deaths / (ifr))
            #sim_recovered = sum(param_values[sim, param_names])
            
            #if (new_recovered > sim_recovered){
            #    difference = new_recovered - sim_recovered
                
            #    # draw from susceptible and put into recovered
            #    probs = param_values[sim, S_names]
            #    draw = rmultinom(1, difference, probs)
            #    param_values[sim, S_names] = param_values[sim, S_names] - draw
            #    stopifnot(all(param_values[sim, S_names] >= 0))
            #    param_values[sim, param_names] = param_values[sim, param_names] + draw
                
            #} else{
            #    difference = sim_recovered - new_recovered
            #    
            #    # draw from recovered and put into susceptible
            #    probs = param_values[sim, param_names]
            #    draw = rmultinom(1, difference, probs)
            #    param_values[sim, param_names] = param_values[sim, param_names] - draw
            #    stopifnot(all(param_values[sim, param_names] >= 0))
            #    param_values[sim, S_names] = param_values[sim, S_names] + draw
            #}
        }
    }
    param_values
}


add_R0_to_output = function(plotdf, R0_vals, parset_id){

foreach(rnum=0:11, .combine='rbind') %do%{
    
    region_to_plot = rnum
    if (rnum==0){
        pop_total = sum(population_list$POPULATION)
    } else{
        pop_total = population_list %>% filter(covid_region == rnum)
        pop_total = sum(pop_total$POPULATION)        
    }

    plotdf %>% 
        filter(Compartment %in% c('S')) %>% 
        mutate(Compartment = case_when((Compartment == 'S') ~ 'Susceptible')) %>%
        filter(as.numeric(Region)==region_to_plot) %>%
        mutate(Cases = Cases / pop_total) -> incidence_plots

        Rt = left_join(incidence_plots, R0_vals, by=c("Date", "Region", "parset")) %>%
            mutate(Rt = Cases * R0) %>% 
            select(-Cases, -Compartment) %>%
            mutate(Cases=Rt, Compartment='Rt') %>%
            select(-Rt, -time, -R0)

    } -> Rt_values

    return(Rt_values)

}