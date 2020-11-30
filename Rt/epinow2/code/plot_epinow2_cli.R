library(readr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())



## Compare outputs from 10-28 (last week) with this week's preliminary outputs

rr <- readr::read_csv('../epinow2_cli_estimates/2020-11-04_EpiNow2v1.2.1_cli_exact/Rt_outputs.csv',
                      col_types = cols(
                        region = col_character(),
                        date = col_date(format = ""),
                        variable = col_character(),
                        strat = col_logical(),
                        type = col_character(),
                        median = col_double(),
                        mean = col_double(),
                        sd = col_double(),
                        lower_95 = col_double(),
                        lower_90 = col_double(),
                        lower_80 = col_double(),
                        upper_80 = col_double(),
                        upper_90 = col_double(),
                        upper_95 = col_double()
                      )) %>%
  select(date, region, median, lower_95, upper_95)

parse_region_names <- function(region){
  ifelse(grepl(region, pattern = '(covidregion)'), gsub(region, pattern = 'covidregion_(\\d\\d?)', replacement = '\\1'), region)
}

rr_last_week <- read_csv('../epinow2_cli_estimates/2020-10-28/epinow2_estimates.csv') %>%
  mutate(region = parse_region_names(geography_modelled)) %>%
  rename(lower_95 = rt_lower,
         median = rt_median,
         upper_95 = rt_upper) %>%
  select(date, region, median, lower_95, upper_95)

bind_rows(list(`2020-10-28` = rr_last_week, `2020-11-04` = rr), .id = 'analysis_date') %>% 
  ggplot() +
  geom_line(aes(x = date, y = median, color = analysis_date))+
  geom_ribbon(aes(x = date, ymin = lower_95, ymax = upper_95, fill = analysis_date), alpha = .5)+
  facet_wrap(.~region)+
  ylab('Rt from CLI data')+
  scale_fill_discrete('')+scale_color_discrete('') ## Remove "kind" header from legend
