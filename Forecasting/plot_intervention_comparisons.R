
get_summary_df = function(input_filenames = input_files, 
                          start_date = as.Date("2020-03-13"),
                          end_date = as.Date("2020-07-31"), 
                          vent_frac = 0.74,
                          save_file = F){
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
  
  df <- data.frame()
  for(i in c(1:length(input_filenames))){
    print(i)
    df_input_sub <- read.csv(input_filenames[i])
    df <- rbind(df,df_input_sub)
  }
  
  df_input <- df

# All infections 
df_infections_summary <- df_input %>%
    mutate(Date = as.Date(Date)) %>% 
    filter(Date < end_date) %>% 
    ungroup() %>% 
    group_by(Time, Compartment) %>%
    filter(Compartment %in% c('E', 'A', 'P', 'IS', 'IM', 'IH', 'IC')) %>%
    group_by(Parset,SimID, Date) %>%
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
  mutate(Date = as.Date(Date)) %>% 
  filter(Date < end_date) %>% 
  ungroup() %>% 
  group_by(Time, Compartment) %>%
  filter(Compartment %in% c('Incidence')) %>%
  group_by(Parset,SimID, Date) %>%
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
  mutate(Date = as.Date(Date)) %>% 
  filter(Date < end_date) %>% 
  ungroup() %>% 
  group_by(Time, Compartment) %>%
  filter(Compartment %in% c('S')) %>%
  group_by(Parset,SimID, Date) %>%
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

# New deaths
df_deaths_summary <- df_input %>%
  mutate(Date = as.Date(Date)) %>% 
  filter(Date < end_date) %>% 
  ungroup() %>% 
  group_by(Time, Compartment) %>%
  filter(Compartment %in% c('nD')) %>%
  group_by(Parset,SimID, Date) %>%
  arrange(Time) %>%
  summarise(Compartment = 'nD',
            Cases = sum(Cases)) %>%
  ungroup() %>% 
  group_by(Date) %>% 
  summarise(
    median = quantile(Cases, 0.5),
    UCI = quantile(Cases, 0.025),
    LCI = quantile(Cases, 0.975)) %>%
  ungroup()%>% 
  mutate(Compartment = "New observed and unobserved deaths")

# Hospitalizations
df_hosp_summary <- df_input %>%
  mutate(Date = as.Date(Date)) %>% 
  filter(Date < end_date) %>% 
  ungroup() %>% 
  group_by(Time, Compartment) %>%
  filter(Compartment %in% c('IH')) %>%
  group_by(Parset,SimID, Date) %>%
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
  mutate(Date = as.Date(Date)) %>% 
  filter(Date < end_date) %>% 
  ungroup() %>% 
  group_by(Time, Compartment) %>%
  filter(Compartment %in% c('IC')) %>%
  group_by(Parset,SimID, Date) %>%
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

# Ventillators in the ICU
df_vent_summary <- df_ICU_summary %>% 
  mutate(median = vent_frac*median, 
        LCI = vent_frac*median,
        UCI = vent_frac*UCI,
        Compartment = "Total Ventillators")

# Recoveries
df_recoveries_summary <- df_input %>%
  mutate(Date = as.Date(Date)) %>% 
  filter(Date < end_date) %>% 
  ungroup() %>% 
  group_by(Time, Compartment) %>%
  filter(Compartment %in% c('R')) %>%
  group_by(Parset,SimID, Date) %>%
  arrange(Time) %>%
  summarise(Compartment = 'R',
            Cases = sum(Cases)) %>%
  ungroup() %>% 
  group_by(Date) %>% 
  summarise(
    median = quantile(Cases, 0.5),
    LCI = quantile(Cases, 0.025),
    UCI = quantile(Cases, 0.975)) %>%
  ungroup()%>% 
  mutate(Compartment = "Recoveries")

summary_df <- rbind(df_inc_summary,
                    df_infections_summary,
                    df_deaths_summary,
                    df_susc_summary,
                    df_recoveries_summary,
                    df_hosp_summary,
                    df_ICU_summary,
                    df_vent_summary)


df_peak <- summary_df %>% filter(Compartment %in% c("New observed and unobserved infections",
                            "Prevalence",
                            "New observed and unobserved deaths",
                            "Total hospital beds occupied",
                            "Total ICU beds occupied",
                            "Total Ventillators")) %>% 
  group_by(Compartment) %>% 
  filter(median == max(median))


return(list(summary = summary_df, 
            peak = df_peak
            ))
}


summary_df <- data.frame()
peak_df <- data.frame()
for( i in c(1:length(intervention_scenarios))){
  intervention = intervention_scenarios[i]
  input_files = paste0(file_path,list.files(pattern = paste0("projections_",intervention), path = file_path))
  summary_overall <- get_summary_df(input_filename = input_files)
  summary_df_sub <- as.data.frame(summary_overall$summary)  %>% mutate(intervention = intervention)
  peak_df_sub <- as.data.frame(summary_overall$peak)  %>% mutate(intervention = intervention)
  summary_df <- rbind(summary_df, summary_df_sub)
  peak_df <- rbind(peak_df, peak_df_sub)
}
# Make plots -------------------------

pop_total <- sum(population1$POPULATION + population2$POPULATION + population3$POPULATION)
df_sum <- as.data.frame(summary_df)
df_sum[df_sum$Compartment %in% c("Recoveries", "Susceptibles","Prevalence"),]$median <- df_sum[df_sum$Compartment %in% c("Recoveries", "Susceptibles","Prevalence"),]$median/pop_total
df_sum[df_sum$Compartment %in% c("Recoveries", "Susceptibles","Prevalence"),]$LCI <- df_sum[df_sum$Compartment %in% c("Recoveries", "Susceptibles","Prevalence"),]$LCI/pop_total
df_sum[df_sum$Compartment %in% c("Recoveries", "Susceptibles","Prevalence"),]$UCI <- df_sum[df_sum$Compartment %in% c("Recoveries", "Susceptibles","Prevalence"),]$UCI/pop_total


p_summary_1 <- ggplot(df_sum %>% filter(Compartment %in% c("New observed and unobserved infections", 
                                                          "New observed and unobserved deaths", 
                                                          "Total hospital beds occupied", 
                                                          "Total ICU beds occupied"
                                                          )), 
                    aes(x = Date, y= median)) + geom_line(aes(color = intervention), linetype = 2) + 
  geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
  #scale_y_continuous(trans = log_trans()) + 
                    # breaks = c(0,100,250,500,1000,5000,10000,100000,500000)) + 
  facet_wrap(~Compartment, scales = "free_y") + 
  ylab("") + 
  #ylab("Prevalent observed and unobserved infections") + 
  theme_classic() +
  theme( axis.text = element_text(size = rel(1.2))) + 
  theme(strip.text.x = element_text(size = 11))

save_plot(paste0("summary_1_",plot_filename), p_summary_1, base_width = 12, base_height = 8) 

p_summary_2 <- ggplot(df_sum %>% filter(Compartment %in% c("Susceptibles",
                                                           "Prevalence")), 
aes(x = Date, y= median)) + geom_line(aes(color = intervention), linetype = 2) + 
  geom_ribbon(aes(ymin = LCI, ymax = UCI, color = intervention, fill = intervention), alpha  = 0.2) + 
  ylim(0,1) + 
  facet_wrap(~Compartment) + 
  ylab("") + 
  #ylab("Prevalent observed and unobserved infections") + 
  theme_classic() +
  theme( axis.text = element_text(size = rel(1.2))) + 
  theme(strip.text.x = element_text(size = 11))

save_plot(paste0("summary_2_",plot_filename), p_summary_1, base_width = 12, base_height = 8) 


