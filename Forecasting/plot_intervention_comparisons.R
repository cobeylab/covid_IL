
# Make population a vector that has the regions concatenated together
population1 = read.csv(population_filename_1)
colnames(population1) <- c("AGE_GROUP_MIN", "POPULATION")

population2 = read.csv(population_filename_2)
colnames(population2) <- c("AGE_GROUP_MIN", "POPULATION")

population3 = read.csv(population_filename_3)
colnames(population3) <- c("AGE_GROUP_MIN", "POPULATION")

# import fraction underreported covariate
fraction_underreported = read.csv(fraction_underreported_file)
row.names(fraction_underreported) = fraction_underreported$time


get_summary_df = function(input_filename = input_file, 
                          start_date = as.Date("2020-03-13"),
                          end_date = as.Date("2020-06-13"), 
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
  input_file = paste0(output_path,list.files(pattern = paste0("plotting_output_",intervention), path = output_path))
  summary_overall <- get_summary_df(input_filename = input_file, regional_aggregation = regional_aggregation)
  summary_df_sub <- as.data.frame(summary_overall$summary)  %>% mutate(intervention = intervention)
  summary_df <- rbind(summary_df, summary_df_sub)
}
# 
summary_df[summary_df$intervention == "100",]$intervention = "Baseline"
summary_df[summary_df$intervention == "80",]$intervention = "20% increase in transmission"
summary_df[summary_df$intervention == "60",]$intervention = "40% increase in transmission"
summary_df[summary_df$intervention == "40",]$intervention = "60% increase in transmission"


save(summary_df, file = paste0(output_path,"summary.rda"))

# Make plots -------------------------
if(regional_aggregation == T){

  pop_total <- sum(population1$POPULATION + population2$POPULATION + population3$POPULATION)
  summary_df_deaths <- summary_df %>% filter(Compartment == "Daily reported deaths")
  df_sum <- as.data.frame(summary_df)
  df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$median <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$median/pop_total
  df_sum[df_sum$Compartment %in% c( "Susceptibles","Prevalence"),]$LCI <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$LCI/pop_total
  df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$UCI <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$UCI/pop_total
  
  
  df_sum_post <- df_sum %>% filter(Date >= today)
  df_sum_pre <- df_sum %>% filter(Date < today) %>% 
    filter(intervention == "Baseline")
  df_sum_new <- rbind(df_sum_pre, df_sum_post) %>% filter(Compartment %in% c("New observed and unobserved infections", 
                                                                             "Prevalence", 
                                                                             "Total hospital beds occupied", 
                                                                             "Total ICU beds occupied"))
  
  
  df_deaths_post <- summary_df_deaths %>% filter(Date >= today)
  df_deaths_pre <- summary_df_deaths %>% filter(Date < today) %>% 
    filter(intervention == "Baseline")
  df_deaths_new <- rbind(df_deaths_pre, df_deaths_post)
  
  ann_text <- data.frame(Compartment = unique(df_sum_new$Compartment),
                           label = c("","","state capacity", "state capacity"),
                           x = c(as.Date("2020-03-15"),as.Date("2020-03-15"), as.Date("2020-04-01"), as.Date("2020-04-01")),
                           y = c(0,0, hospital_capacity + hospital_capacity/ICU_capacity*200, ICU_capacity + 200 )
                                 )
   
  
  df_sum_new$capacity = NA
  df_sum_new[df_sum_new$Compartment == "Total hospital beds occupied",]$capacity = hospital_capacity
  df_sum_new[df_sum_new$Compartment == "Total ICU beds occupied",]$capacity = ICU_capacity
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
    geom_text(data = ann_text,
              mapping = aes(x = x, y = y, label = label))
  
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
    theme(strip.text.x = element_text(size = 10))
  
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
    theme(strip.text.x = element_text(size = 10))
  
  p_output_indefinite <- plot_grid(title, 
                          p_summary_indefinite,
                          ncol = 1,
                          rel_heights = c(0.1,1))
  
  
  save_plot(paste0("summary_indefinite_statewide_",plot_filename), p_output_indefinite, base_width = 12, base_height = 8) 
  
  death_data <- read.csv(IDPH_death_data_filename) %>% mutate(Date = as.Date(date, format = "%m/%d/%Y")) %>% 
    select(Date, new_deaths) %>% 
    mutate(Compartment = "Daily reported deaths",
           intervention = "Reported deaths, IDPH")
  
  p_deaths <- ggplot(df_deaths_new %>% mutate(Compartment = "Daily reported deaths"), aes(x = Date, y = median)) + geom_line(aes(color = intervention), linetype = 2) + 
    geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
    geom_point(data = death_data, aes(x = Date, y = new_deaths)) + 
    facet_wrap(~Compartment, scales = "free_y") + 
    ylab("") + 
    #ylab("Prevalent observed and unobserved infections") + 
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10)) + 
    theme(strip.text.x = element_text(size = 10))
  
  
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
  
  
  p_output_deaths <- plot_grid(title2, 
                              p_deaths,
                              ncol = 1,
                              rel_heights = c(0.1,1))
  
  save_plot("death_summary_statewide_outputs.png",p_output_deaths, base_width = 8, base_height = 4 )
}


if(regional_aggregation == F){
  regions = unique(summary_df$Region)
  population_vec <- c(sum(population1$POPULATION), sum(population2$POPULATION), sum(population3$POPULATION))
  for(region in regions){
    pop_total <- population_vec[region]
    df_sum <- as.data.frame(summary_df %>% filter(Region == region))
    summary_df_deaths <- summary_df %>% filter(Compartment == "Daily reported deaths" & Region == region)
   
    
    df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$median <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$median/pop_total
    df_sum[df_sum$Compartment %in% c( "Susceptibles","Prevalence"),]$LCI <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$LCI/pop_total
    df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$UCI <- df_sum[df_sum$Compartment %in% c("Susceptibles","Prevalence"),]$UCI/pop_total
    
    
    df_sum_post <- df_sum %>% filter(Date >= today)
    df_sum_pre <- df_sum %>% filter(Date < today) %>% 
      filter(intervention == "Baseline")
    df_sum_new <- rbind(df_sum_pre, df_sum_post) %>% filter(Compartment %in% c("New observed and unobserved infections", 
                                                                               "Prevalence", 
                                                                               "Total hospital beds occupied", 
                                                                               "Total ICU beds occupied"))
    
    
    df_deaths_post <- summary_df_deaths %>% filter(Date >= today)
    df_deaths_pre <- summary_df_deaths %>% filter(Date < today) %>% 
      filter(intervention == "Baseline")
    df_deaths_new <- rbind(df_deaths_pre, df_deaths_post)
    
    ann_text <- data.frame(Compartment = unique(df_sum_new$Compartment),
                           label = c("","","state capacity", "state capacity"),
                           x = c(as.Date("2020-03-15"),as.Date("2020-03-15"), as.Date("2020-04-01"), as.Date("2020-04-01")),
                           y = c(0,0, hospital_capacity + hospital_capacity/ICU_capacity*200, ICU_capacity + 200 )
    )
    
    
    df_sum_new$capacity = NA
    df_sum_new[df_sum_new$Compartment == "Total hospital beds occupied",]$capacity = hospital_capacity
    df_sum_new[df_sum_new$Compartment == "Total ICU beds occupied",]$capacity = ICU_capacity
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
      geom_text(data = ann_text,
                mapping = aes(x = x, y = y, label = label))
    
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
      theme(strip.text.x = element_text(size = 10))
    
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
      theme(strip.text.x = element_text(size = 10))
    
    p_output_indefinite <- plot_grid(title, 
                                     p_summary_indefinite,
                                     ncol = 1,
                                     rel_heights = c(0.1,1))
    
    
    save_plot(paste0("region_", region, "_summary_indefinite_",plot_filename), p_output_indefinite, base_width = 12, base_height = 8) 
    
    p_deaths <- ggplot(df_deaths_new %>% mutate(Compartment = "Daily reported deaths"), aes(x = Date, y = median)) + geom_line(aes(color = intervention), linetype = 2) + 
      geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
      facet_wrap(~Compartment, scales = "free_y") + 
      ylab("") + 
      #ylab("Prevalent observed and unobserved infections") + 
      theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10)) + 
      theme(strip.text.x = element_text(size = 10))
    
    
    title2 <- ggdraw() + 
      draw_label(
        paste0("Simulated reported deaths (median and 95% CI), data through ", this_date),
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



                      