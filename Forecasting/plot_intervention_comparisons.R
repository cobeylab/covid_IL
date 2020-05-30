
# Make population a vector that has the regions concatenated together
population1 = read.csv(population_filename_1)
colnames(population1) <- c("AGE_GROUP_MIN", "POPULATION")

population2 = read.csv(population_filename_2)
colnames(population2) <- c("AGE_GROUP_MIN", "POPULATION")

population3 = read.csv(population_filename_3)
colnames(population3) <- c("AGE_GROUP_MIN", "POPULATION")

population4 = read.csv(population_filename_4)
colnames(population4) <- c("AGE_GROUP_MIN", "POPULATION")

# import fraction underreported covariate
fraction_underreported = read.csv(fraction_underreported_file)
row.names(fraction_underreported) = fraction_underreported$time


get_summary_df = function(input_filename = input_file, 
                          start_date = as.Date("2020-03-13"),
                          end_date = as.Date("2020-10-01"), 
                          regional_aggregation,
                          vent_frac = 0.74){
  ## Load packages and specify files 
  require(dplyr)
  require(pomp)
  require(ggplot2)
  require(magrittr)
  require(lubridate)
  require(MASS)
  require(reshape2)
  require(dplyr)
  select <- dplyr::select
  rename <- dplyr::rename
  summarize <- dplyr::summarise
  contains <- dplyr::contains

  print(input_filename)
  df <-  read.csv(input_filename)

  
  df_input <- df  %>% 
    mutate(Date = as.Date(Date)) %>% 
    filter(Date >= start_date & Date <= end_date) 

  
if(regional_aggregation == T){
  # All infections 
  df_infections_summary <- df_input %>%
      ungroup() %>% 
      group_by(Time, Compartment) %>%
      filter(Compartment %in% c('E', 'A', 'P', 'IS', 'IM', 'IH', 'IC')) %>%
      group_by(parset,SimID, Date) %>%
      arrange(Time) %>%
      summarise(Compartment = 'I',
                Cases = sum(Cases)) %>%
      ungroup() %>% 
      group_by(Date) %>% 
      summarise(
        median = quantile(Cases, 0.5),
        LCI = quantile(Cases, 0.025),
        UCI = quantile(Cases, 0.975)) %>%
    ungroup() %>% 
    mutate(Compartment = "Prevalence")
  
  ## New infections
  df_inc_summary <- df_input %>%
    ungroup() %>% 
    group_by(Time, Compartment) %>%
    filter(Compartment %in% c('Incidence')) %>%
    group_by(parset,SimID, Date) %>%
    arrange(Time) %>%
    summarise(Compartment = 'I',
              Cases = sum(Cases)) %>%
    ungroup() %>% 
    group_by(Date) %>% 
    summarise(
      median = quantile(Cases, 0.5),
      LCI = quantile(Cases, 0.025),
      UCI = quantile(Cases, 0.975)) %>%
    ungroup() %>% 
    mutate(Compartment = "New observed and unobserved infections")
  
  ## Susceptibles
  df_susc_summary <- df_input %>%
    ungroup() %>% 
    group_by(Time, Compartment) %>%
    filter(Compartment %in% c('S')) %>%
    group_by(parset,SimID, Date) %>%
    arrange(Time) %>%
    summarise(Compartment = 'S',
              Cases = sum(Cases)) %>%
    ungroup() %>% 
    group_by(Date) %>% 
    summarise(
      median = quantile(Cases, 0.5),
      LCI = quantile(Cases, 0.025),
      UCI = quantile(Cases, 0.975)) %>%
    ungroup()%>% 
    mutate(Compartment = "Susceptibles")
  
  
  # Hospitalizations
  df_hosp_summary <- df_input %>%
    ungroup() %>% 
    group_by(Time, Compartment) %>%
    filter(Compartment %in% c('IH','IC')) %>%
    group_by(parset,SimID, Date) %>%
    arrange(Time) %>%
    summarise(Compartment = 'IH',
              Cases = sum(Cases)) %>%
    ungroup() %>% 
    group_by(Date) %>% 
    summarise(
      median = quantile(Cases, 0.5),
      LCI = quantile(Cases, 0.025),
      UCI = quantile(Cases, 0.975)) %>%
    ungroup() %>% 
    mutate(Compartment = "Total hospital beds occupied")
  
  # ICU beds
  df_ICU_summary <- df_input %>%
    ungroup() %>% 
    group_by(Time, Compartment) %>%
    filter(Compartment %in% c('IC')) %>%
    group_by(parset,SimID, Date) %>%
    arrange(Time) %>%
    summarise(Compartment = 'IC',
              Cases = sum(Cases)) %>%
    ungroup() %>% 
    group_by(Date) %>% 
    summarise(
      median = quantile(Cases, 0.5),
      LCI = quantile(Cases, 0.025),
      UCI = quantile(Cases, 0.975)) %>%
    ungroup()%>% 
    mutate(Compartment = "Total ICU beds occupied")
  
  # ICU beds
  df_reported_deaths_summary <- df_input %>%
    ungroup() %>% 
    group_by(Time, Compartment) %>%
    filter(Compartment %in% c('Reported hd', 'Reported nhd')) %>%
    group_by(parset,SimID, Date) %>%
    arrange(Time) %>%
    summarise(Compartment = 'RD',
              Cases = sum(Cases)) %>%
    ungroup() %>% 
    group_by(Date) %>% 
    summarise(
      median = quantile(Cases, 0.5),
      LCI = quantile(Cases, 0.025),
      UCI = quantile(Cases, 0.975)) %>%
    ungroup()%>% 
    mutate(Compartment = "Daily reported deaths")
} 

  
  if(regional_aggregation == F){
    # All infections 
    df_infections_summary <- df_input %>%
      ungroup() %>% 
      group_by(Time, Compartment, Region) %>%
      filter(Compartment %in% c('E', 'A', 'P', 'IS', 'IM', 'IH', 'IC')) %>%
      group_by(parset,SimID, Date, Region) %>%
      arrange(Time) %>%
      summarise(Compartment = 'I',
                Cases = sum(Cases)) %>%
      ungroup() %>% 
      group_by(Date,Region) %>% 
      summarise(
        median = quantile(Cases, 0.5),
        LCI = quantile(Cases, 0.025),
        UCI = quantile(Cases, 0.975)) %>%
      ungroup() %>% 
      mutate(Compartment = "Prevalence")
    
    ## New infections
    df_inc_summary <- df_input %>%
      ungroup() %>% 
      group_by(Time, Compartment, Region) %>%
      filter(Compartment %in% c('Incidence')) %>%
      group_by(parset,SimID, Date, Region) %>%
      arrange(Time) %>%
      summarise(Compartment = 'I',
                Cases = sum(Cases)) %>%
      ungroup() %>% 
      group_by(Date, Region) %>% 
      summarise(
        median = quantile(Cases, 0.5),
        LCI = quantile(Cases, 0.025),
        UCI = quantile(Cases, 0.975)) %>%
      ungroup() %>% 
      mutate(Compartment = "New observed and unobserved infections")
    
    ## Susceptibles
    df_susc_summary <- df_input %>%
      ungroup() %>% 
      group_by(Time, Compartment, Region) %>%
      filter(Compartment %in% c('S')) %>%
      group_by(parset,SimID, Date, Region) %>%
      arrange(Time) %>%
      summarise(Compartment = 'S',
                Cases = sum(Cases)) %>%
      ungroup() %>% 
      group_by(Date, Region) %>% 
      summarise(
        median = quantile(Cases, 0.5),
        LCI = quantile(Cases, 0.025),
        UCI = quantile(Cases, 0.975)) %>%
      ungroup()%>% 
      mutate(Compartment = "Susceptibles")
    
    
    # Hospitalizations
    df_hosp_summary <- df_input %>%
      ungroup() %>% 
      group_by(Time, Compartment, Region) %>%
      filter(Compartment %in% c('IH','IC')) %>%
      group_by(parset, SimID, Date, Region) %>%
      arrange(Time) %>%
      summarise(Compartment = 'IH',
                Cases = sum(Cases)) %>%
      ungroup() %>% 
      group_by(Date, Region) %>% 
      summarise(
        median = quantile(Cases, 0.5),
        LCI = quantile(Cases, 0.025),
        UCI = quantile(Cases, 0.975)) %>%
      ungroup() %>% 
      mutate(Compartment = "Total hospital beds occupied")
    
    # ICU beds
    df_ICU_summary <- df_input %>%
      ungroup() %>% 
      group_by(Time, Compartment, Region) %>%
      filter(Compartment %in% c('IC')) %>%
      group_by(parset,SimID, Date, Region) %>%
      arrange(Time) %>%
      summarise(Compartment = 'IC',
                Cases = sum(Cases)) %>%
      ungroup() %>% 
      group_by(Date, Region) %>% 
      summarise(
        median = quantile(Cases, 0.5),
        LCI = quantile(Cases, 0.025),
        UCI = quantile(Cases, 0.975)) %>%
      ungroup()%>% 
      mutate(Compartment = "Total ICU beds occupied")
    
    # ICU beds
    df_reported_deaths_summary <- df_input %>%
      ungroup() %>% 
      group_by(Time, Compartment, Region) %>%
      filter(Compartment %in% c('Reported hd', 'Reported nhd')) %>%
      group_by(parset,SimID, Date, Region) %>%
      arrange(Time) %>%
      summarise(Compartment = 'RD',
                Cases = sum(Cases)) %>%
      ungroup() %>% 
      group_by(Date, Region) %>% 
      summarise(
        median = quantile(Cases, 0.5),
        LCI = quantile(Cases, 0.025),
        UCI = quantile(Cases, 0.975)) %>%
      ungroup()%>% 
      mutate(Compartment = "Daily reported deaths")
  }  
  summary_df <- rbind(df_inc_summary,
                      df_infections_summary,
                      df_susc_summary,
                      df_hosp_summary,
                      df_ICU_summary,
                      df_reported_deaths_summary)


return(list(summary = summary_df
            ))
}

summary_df <- data.frame()
for( i in c(1:length(intervention_scenarios))){
  intervention = intervention_scenarios[i]

  if(regional_aggregation){
      input_file = paste0(output_path,sprintf("plotting_output_statewide_%s.csv", intervention))
  } else{
      input_file = paste0(output_path,sprintf("plotting_output_%s.csv", intervention))
  }

  print(input_file)
  summary_overall <- get_summary_df(input_filename = input_file, regional_aggregation = regional_aggregation)
  summary_df_sub <- as.data.frame(summary_overall$summary)  %>% mutate(intervention = intervention)
  summary_df <- rbind(summary_df, summary_df_sub)
}
# 
summary_df[summary_df$intervention == "100",]$intervention = "Baseline"
summary_df[summary_df$intervention == "10",]$intervention = "10% of the way to pre-SIP"
summary_df[summary_df$intervention == "30",]$intervention = "30% of the way to pre-SIP"
#summary_df[summary_df$intervention == "80",]$intervention = "20% increase in transmission"
#summary_df[summary_df$intervention == "60",]$intervention = "40% increase in transmission"
#summary_df[summary_df$intervention == "40",]$intervention = "60% increase in transmission"


save(summary_df, file = paste0(output_path,"summary.rda"))

# Make plots -------------------------
if(regional_aggregation == T){

    
    hosp_cap = reg_hosp_cap %>% 
         mutate(date = as.Date(date)) %>% 
        filter(date == max(date)) %>% 
        group_by(date) %>%
        summarize(capacity_non_icu=sum(capacity_non_icu),
                  capacity_adult_icu=sum(capacity_adult_icu))
    hospital_capacity = as.numeric(hosp_cap$capacity_non_icu + hosp_cap$capacity_adult_icu)
    ICU_capacity = as.numeric(hosp_cap$capacity_adult_icu)

  pop_total <- sum(population1$POPULATION + population2$POPULATION + population3$POPULATION + population4$POPULATION)
  summary_df_deaths <- summary_df %>% filter(Compartment == "Daily reported deaths")
  df_sum <- as.data.frame(summary_df)
  df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$median <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$median/pop_total
  df_sum[df_sum$Compartment %in% c( "Susceptibles","Prevalence"),]$LCI <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$LCI/pop_total
  df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$UCI <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$UCI/pop_total
  
  
  df_sum_post <- df_sum %>% filter(Date >= intervention_lift_date) #today)
  df_sum_pre <- df_sum %>% filter(Date < intervention_lift_date) %>% #today) %>% 
    filter(intervention == "Baseline")
  df_sum_new <- rbind(df_sum_pre, df_sum_post) %>% filter(Compartment %in% c("New observed and unobserved infections", 
                                                                             "Prevalence", 
                                                                             "Total hospital beds occupied", 
                                                                             "Total ICU beds occupied"))
  
  
  df_deaths_post <- summary_df_deaths %>% filter(Date >= intervention_lift_date)
  df_deaths_pre <- summary_df_deaths %>% filter(Date < intervention_lift_date) %>% 
    filter(intervention == "Baseline")
  df_deaths_new <- rbind(df_deaths_pre, df_deaths_post)
  
    df_sum_new$capacity=NA
    df_sum_new[df_sum_new$Compartment == "Total hospital beds occupied",]$capacity = hospital_capacity
    df_sum_new[df_sum_new$Compartment == "Total ICU beds occupied",]$capacity = ICU_capacity
    print('madeit')

  limits = df_sum_new %>% group_by(Compartment) %>% summarize(ma = max(UCI, capacity, na.rm=TRUE), mi = min(LCI))
  midpoints = (limits$ma - limits$mi)/2
  
  death_lim = df_deaths_new %>%  summarize(ma=max(UCI), mi=min(LCI))
  death_mid = (death_lim$ma - death_lim$mi)/2
  
  
  ann_text <- data.frame(Compartment = unique(df_sum_new$Compartment),
                           label = c("","","state capacity", "state capacity"),
                           x = c(as.Date("2020-03-15"),as.Date("2020-03-15"), as.Date("2020-04-01"), as.Date("2020-04-01")),
                           y = c(0,0, hospital_capacity + hospital_capacity/ICU_capacity*200, ICU_capacity + 700 )
                                 )
  ann_text_lift = data.frame(Compartment = unique(df_sum_new$Compartment),
                             label = "relax SIP",
                             x = intervention_lift_date - 3,
                             y = midpoints
  )   
  
  ann_text_forecast = data.frame(Compartment = unique(df_sum_new$Compartment),
                             label = "forecast start",
                             x = today - 3,
                             y = midpoints
  )   
  
  ann_text_deaths_lift <- data.frame(Compartment = "Daily reported deaths",
                           label = "relax SIP",
                           x = intervention_lift_date - 5,
                           y = death_mid + 20)

  ann_text_deaths_forecast <- data.frame(Compartment = "Daily reported deaths",
                           label = "forecast start",
                           x = today - 5,
                           y = death_mid)  

  p_summary_1 <- ggplot(df_sum_new,aes(x = Date, y= median)) + 
    geom_line(aes(color = intervention), linetype = 2) + 
    geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
    geom_hline(aes(yintercept = capacity), linetype = 2) + 
    facet_wrap(~Compartment, scales = "free_y") + 
    ylab("") + 
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10)) + 
    theme(strip.text.x = element_text(size = 10)) +
    labs( color='Change in transmission', fill='Change in transmission') +
    geom_text(data = ann_text,
              mapping = aes(x = x, y = y, label = label)) +
    scale_x_date(date_labels = "%b %d") +
    geom_vline(xintercept=today) + 
      geom_text(data = ann_text_lift,
                mapping = aes(x = x, y = y, label = label),
                angle=90) +
      geom_text(data = ann_text_forecast,
                mapping = aes(x = x, y = y, label = label),
                angle=90) +
      geom_vline(xintercept=today) +
      geom_vline(xintercept=intervention_lift_date)
  
  title <- ggdraw() + 
    draw_label(
      paste0("Simulated projections (median and 95% CI) for models fitted to data through ", this_date),
      fontface = 'bold',
      size = 12,
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  p_output_1 <- plot_grid(title, 
                       p_summary_1,
                       ncol = 1,
                       rel_heights = c(0.1,1))
  
  
  save_plot(paste0("summary_1_statewide_",plot_filename), p_output_1, base_width = 12, base_height = 8) 
  
  p_summary_2 <- ggplot(df_sum %>% filter(Compartment %in% c("Susceptibles",
                                                             "Prevalence")), 
  aes(x = Date, y= median)) + geom_line(aes(color = intervention), linetype = 2) + 
    geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
    ylim(0,1) + 
    facet_wrap(~Compartment) + 
    ylab("") + 
    #ylab("Prevalent observed and unobserved infections") + 
    theme_classic() +
    theme( axis.text = element_text(size = rel(1.1))) + 
    theme(strip.text.x = element_text(size = 10)) +
    scale_x_date(date_labels = "%b %d")
  
  p_output_2 <- plot_grid(title, 
                          p_summary_2,
                          ncol = 1,
                          rel_heights = c(0.1,1))
  
  save_plot(paste0("summary_2_statewide_",plot_filename), p_output_2, base_width = 12, base_height = 4) 
  
  
  p_summary_indefinite <- ggplot(df_sum %>% filter(intervention == "Baseline" & Compartment %in% c("New observed and unobserved infections", 
                                                                                                     "Prevalence",
                                                                      "Total hospital beds occupied", 
                                                                      "Total ICU beds occupied")), 
  aes(x = Date, y= median)) + geom_line(aes(color = intervention), linetype = 2) + 
    geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
    #scale_y_continuous(trans = log_trans()) + 
    # breaks = c(0,100,250,500,1000,5000,10000,100000,500000)) + 
    facet_wrap(~Compartment, scales = "free_y") + 
    ylab("") + 
    #ylab("Prevalent observed and unobserved infections") + 
    theme_classic() +
    theme( axis.text = element_text(size = rel(1.2))) + 
    theme(strip.text.x = element_text(size = 10)) +
    scale_x_date(date_labels = "%b %d")
  
  p_output_indefinite <- plot_grid(title, 
                          p_summary_indefinite,
                          ncol = 1,
                          rel_heights = c(0.1,1))
  
  
  save_plot(paste0("summary_indefinite_statewide_",plot_filename), p_output_indefinite, base_width = 12, base_height = 8) 
  
  death_data <- read.csv(IDPH_death_data_filename) %>% mutate(Date = as.Date(date)) %>% 
    select(Date, new_deaths) %>% 
    mutate(Compartment = "Daily reported deaths",
           intervention = "Reported deaths, IDPH") %>%
      filter(Date <= today)
  
 # death_data <- read.csv('civis_data/total_deaths.csv') %>% mutate(Date = as.Date(date)) %>% 
 #     group_by(Date) %>%
 #     summarize(new_deaths = sum(new_deaths)) %>%
 #     ungroup() %>%
 #     mutate(Compartment = "Daily reported deaths",
 #            intervention = "Reported deaths, IDPH") %>%
 #     filter(Date <= today)
  
  p_deaths <- ggplot(df_deaths_new %>% mutate(Compartment = "Daily reported deaths"), aes(x = Date, y = median)) + geom_line(aes(color = intervention), linetype = 2) + 
    geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
    geom_point(data = death_data, aes(x = Date, y = new_deaths)) + 
    #facet_wrap(~Compartment, scales = "free_y") + 
    ylab("") + 
    #ylab("Prevalent observed and unobserved infections") + 
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10)) + 
    theme(strip.text.x = element_text(size = 10)) +
      labs( color='Change in transmission', fill = 'Change in transmission') +
    scale_x_date(date_labels = "%b %d") +
    geom_text(data = ann_text_deaths_lift,
                mapping = aes(x = x, y = y, label = label),
                angle=90) +
    #geom_text(data = ann_text_deaths_forecast,
    #            mapping = aes(x = x, y = y, label = label),
    #            angle=90) +
    #geom_vline(xintercept=today) +
    geom_vline(xintercept=intervention_lift_date) +
    ylab('Daily reported deaths')

  
  
  title2 <- ggdraw() + 
    draw_label(
      paste0("Simulated reported deaths (median and 95% CI) and reported deaths from IDPH (points), data through ", this_date),
      fontface = 'bold',
      size = 10,
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  
  p_output_deaths <- plot_grid(#title2, 
                              p_deaths,
                              ncol = 1,
                              rel_heights = c(0.1,1))
  
  save_plot("death_summary_statewide_outputs.png",p_output_deaths, base_width = 8, base_height = 4 )
}


if(regional_aggregation == F){
    
    all_death_data <- read.csv(regional_idph_file) %>% 
        mutate(Date = as.Date(test_date)) %>% 
        select(Date, inc_deaths, region) %>% 
        mutate(Compartment = "Daily reported deaths",
               intervention = "Reported deaths, IDPH") %>%
        filter(inc_deaths >= 0, Date <= today)
    
    
  regions = unique(summary_df$Region)
  population_vec <- c(sum(population1$POPULATION), sum(population2$POPULATION), sum(population3$POPULATION), sum(population4$POPULATION))
  region_names = c('north-central', 'central', 'northeast', 'southern')
  for(region in regions){

    region_name = region_names[region]
    death_data = all_death_data %>% filter(region==region_name)
    pop_total <- population_vec[region]
    df_sum <- as.data.frame(summary_df %>% filter(Region == region))
    summary_df_deaths <- summary_df %>% filter(Compartment == "Daily reported deaths" & Region == region)
   
    
    df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$median <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$median/pop_total
    df_sum[df_sum$Compartment %in% c( "Susceptibles","Prevalence"),]$LCI <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$LCI/pop_total
    df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$UCI <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$UCI/pop_total
    
    
    hosp_cap = reg_hosp_cap %>% filter(idph_cluster == region_name) %>% mutate(date = as.Date(date)) %>% filter(date == max(date))
    
    df_sum_post <- df_sum %>% filter(Date >= today)
    df_sum_pre <- df_sum %>% filter(Date < today) %>% 
      filter(intervention == "Baseline")
    df_sum_new <- rbind(df_sum_pre, df_sum_post) %>% filter(Compartment %in% c("New observed and unobserved infections", 
                                                                               "Prevalence", 
                                                                               "Total hospital beds occupied", 
                                                                               "Total ICU beds occupied"))
    

    
    df_deaths_post <- summary_df_deaths %>% filter(Date >= intervention_lift_date)
    df_deaths_pre <- summary_df_deaths %>% filter(Date < intervention_lift_date) %>% 
      filter(intervention == "Baseline")
    df_deaths_new <- rbind(df_deaths_pre, df_deaths_post)
    
    hospital_capacity = as.numeric(hosp_cap$capacity_non_icu + hosp_cap$capacity_adult_icu)
    ICU_capacity = as.numeric(hosp_cap$capacity_adult_icu)
    
    df_sum_new$capacity = NA
    df_sum_new[df_sum_new$Compartment == "Total hospital beds occupied",]$capacity = hospital_capacity + ICU_capacity
    df_sum_new[df_sum_new$Compartment == "Total ICU beds occupied",]$capacity = ICU_capacity
    
    limits = df_sum_new %>% group_by(Compartment) %>% summarize(ma = max(UCI, capacity, na.rm=TRUE), mi = min(LCI))
    midpoints = (limits$ma - limits$mi)/2
    
    death_lim = df_deaths_new %>%  summarize(ma=max(UCI), mi=min(LCI))
    death_mid = (death_lim$ma - death_lim$mi)/2
    
    ann_text <- data.frame(Compartment = unique(df_sum_new$Compartment),
                           label = c(NA,NA,"approximate\nregional capacity", "approximate\nregional capacity"),
                           x = c(as.Date("2020-03-15"),as.Date("2020-03-15"), as.Date("2020-04-05"), as.Date("2020-04-05")),
                           y = c(0,0, hospital_capacity + ICU_capacity, ICU_capacity)
    )
    
    ann_text_lift = data.frame(Compartment = unique(df_sum_new$Compartment),
                               label = "relax SIP",
                               x = intervention_lift_date - 3,
                               y = midpoints
    )   
    
    ann_text_forecast = data.frame(Compartment = unique(df_sum_new$Compartment),
                                   label = "forecast start",
                                   x = today - 3,
                                   y = midpoints
    )   
    
    ann_text_deaths_lift <- data.frame(Compartment = "Daily reported deaths",
                                       label = "relax SIP",
                                       x = intervention_lift_date - 3,
                                       y = death_mid)
    
    ann_text_deaths_forecast <- data.frame(Compartment = "Daily reported deaths",
                                           label = "forecast start",
                                           x = today - 3,
                                           y = death_mid) 
    

    p_summary_1 <- ggplot(df_sum_new,aes(x = Date, y= median)) + 
      geom_line(aes(color = intervention), linetype = 2) + 
      geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) +
      geom_hline(aes(yintercept=capacity), linetype = 2) +
        geom_label(data = ann_text,
                  mapping = aes(x = x, y = y, label = label), label.size=0) +
        geom_text(data = ann_text_lift,
                  mapping = aes(x = x, y = y, label = label),
                  angle=90) +
        geom_text(data = ann_text_forecast,
                  mapping = aes(x = x, y = y, label = label),
                  angle=90) +
        geom_vline(xintercept=today) +
        geom_vline(xintercept=intervention_lift_date) +
      facet_wrap(~Compartment, scales = "free_y") + 
      ylab("") + 
      theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10)) + 
      theme(strip.text.x = element_text(size = 10)) + 
      scale_x_date(date_labels = "%b %d") +
        labs( color='Change in transmission', fill='Change in transmission')
    
    title <- ggdraw() + 
      draw_label(
        sprintf("Simulated projections (median and 95%% CI) for %s region; models fitted to data through %s", region_name, this_date),
        fontface = 'bold',
        size = 12,
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
    
    p_output_1 <- plot_grid(title, 
                            p_summary_1,
                            ncol = 1,
                            rel_heights = c(0.1,1))
    
    
    save_plot(paste0("region_", region, "_summary_1__",plot_filename), p_output_1, base_width = 12, base_height = 8) 
    
    p_summary_2 <- ggplot(df_sum %>% filter(Compartment %in% c("Susceptibles",
                                                               "Prevalence")), 
                          aes(x = Date, y= median)) + geom_line(aes(color = intervention), linetype = 2) + 
      geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
      ylim(0,1) + 
      facet_wrap(~Compartment) + 
      ylab("") + 
      #ylab("Prevalent observed and unobserved infections") + 
      theme_classic() +
      theme( axis.text = element_text(size = rel(1.1))) + 
      theme(strip.text.x = element_text(size = 10)) +
      scale_x_date(date_labels = "%b %d")
    
    p_output_2 <- plot_grid(title, 
                            p_summary_2,
                            ncol = 1,
                            rel_heights = c(0.1,1))
    
    save_plot(paste0("region_", region, "_summary_2_",plot_filename), p_output_2, base_width = 12, base_height = 4) 
    
    
    p_summary_indefinite <- ggplot(df_sum %>% filter(intervention == "Baseline" & Compartment %in% c("New observed and unobserved infections", 
                                                                                                     "Prevalence",
                                                                                                     "Total hospital beds occupied", 
                                                                                                     "Total ICU beds occupied")), 
                                   aes(x = Date, y= median)) + geom_line(aes(color = intervention), linetype = 2) + 
      geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
      #scale_y_continuous(trans = log_trans()) + 
      # breaks = c(0,100,250,500,1000,5000,10000,100000,500000)) + 
      facet_wrap(~Compartment, scales = "free_y") + 
      ylab("") + 
      #ylab("Prevalent observed and unobserved infections") + 
      theme_classic() +
      theme( axis.text = element_text(size = rel(1.2))) + 
      theme(strip.text.x = element_text(size = 10)) +
      scale_x_date(date_labels = "%b %d")
    
    p_output_indefinite <- plot_grid(title, 
                                     p_summary_indefinite,
                                     ncol = 1,
                                     rel_heights = c(0.1,1))
    
    
    save_plot(paste0("region_", region, "_summary_indefinite_",plot_filename), p_output_indefinite, base_width = 12, base_height = 8) 
    
    p_deaths <- ggplot(df_deaths_new %>% mutate(Compartment = "Daily reported deaths"), aes(x = Date, y = median)) + geom_line(aes(color = intervention), linetype = 2) + 
      geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
      geom_point(data=death_data, aes(x=Date, y=inc_deaths)) +
      facet_wrap(~Compartment, scales = "free_y") + 
      ylab("") + 
      #ylab("Prevalent observed and unobserved infections") + 
      theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10)) + 
      theme(strip.text.x = element_text(size = 10)) +
      scale_x_date(date_labels = "%b %d") +
        geom_text(data = ann_text_deaths_lift,
                  mapping = aes(x = x, y = y, label = label),
                  angle=90) +
        geom_text(data = ann_text_deaths_forecast,
                  mapping = aes(x = x, y = y, label = label),
                  angle=90) +
        geom_vline(xintercept=today) +
        geom_vline(xintercept=intervention_lift_date) +
        labs( color='Change in transmission', fill='Change in transmission')
    
    
    
    title2 <- ggdraw() + 
      draw_label(
          sprintf("Simulated reported deaths (median and 95%% CI) for %s region; models fitted to data through %s", region_name, this_date),
        fontface = 'bold',
        size = 10,
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
    
    
    p_output_deaths <- plot_grid(title2, 
                                 p_deaths,
                                 ncol = 1,
                                 rel_heights = c(0.1,1))
    
    save_plot(paste0("region_",region, "_death_summary_outputs.png"),p_output_deaths, base_width = 8, base_height = 4 )
  }
}



                      