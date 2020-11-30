min_0 <- function(xx){ifelse(xx<0, 0, xx)}

load_cli <- function(){
  ## Load the public linelist data by restore region ------------
  ## More up to date, public linelist
  region_cli = read.csv('../internal_data/cli_admissions_latest.csv') %>%
    mutate(date = as.Date(date)) %>%
    rename(
           nadmit = cli) %>%
    mutate(
           region = as.character(covid_region)) %>%
    mutate(nadmit = ifelse(is.na(nadmit), 0, nadmit)) %>%
    select(-covid_region)

  statewide_cli = region_cli %>%
  group_by(date) %>%
  summarise(region = 'illinois',
            nadmit = sum(nadmit)) %>%
  ungroup()

  bind_rows(region_cli, statewide_cli) %>%
    group_by(region) %>%
    arrange(date) %>%
    mutate(smoothed = smooth.spline(nadmit, spar = .5)$y %>% min_0,
           avg_7d = zoo::rollmean(nadmit, k = 7, fill = c(mean(nadmit[1:7], na.rm = T), 
                                                             NA, 
                                                             mean(nadmit[length(nadmit)-(0:6)], na.rm = T)))
    )%>%
    ungroup() 
}
