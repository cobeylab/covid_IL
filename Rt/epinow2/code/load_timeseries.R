min_0 <- function(xx){ifelse(xx<0, 0, xx)}

load_idph_public_cases_restore_region <- function(){
## Load the public linelist data by restore region ------------
## More up to date, public linelist
raw_dat <- read_csv('../data/idph_public_restore_region.csv') %>%
  rename(region = restore_region) %>%
  filter(region != 'unknown') %>%
  group_by(date, region) %>%
  summarise(new_cases = sum(new_cases)) %>%
  mutate(new_cases = ifelse(is.na(new_cases), 0, new_cases),
         region = toupper(region),
         region = ifelse(region=='NORTH-CENTRAL', 'NORTHCENTRAL', region)) %>%
  ungroup() 
## Get all IL data
overall_dat <- raw_dat %>%
  filter(region != 'CHICAGO') %>%
  group_by(date) %>%
  summarise(region = 'IL_Overall',
            new_cases = sum(new_cases))
## Add IL overall data as its own region
bind_rows(raw_dat, overall_dat) %>%
  ungroup() %>% group_by(region) %>%
  arrange(date) %>%
  mutate(smoothed = smooth.spline(new_cases, spar = .6)$y %>% min_0,
         avg_7d = zoo::rollmean(new_cases, k = 7, fill = c(mean(new_cases[1:7], na.rm = T), 
                                                           NA, 
                                                           mean(new_cases[length(new_cases)-(0:6)], na.rm = T)))
  )%>%
  filter(date <= max(date))%>%
  filter(date >= lubridate::as_date('2020-04-01')) %>%
  mutate(new_cases = ifelse(new_cases<0, 0, new_cases)) %>%
  ungroup()
}





load_idph_public_cases_covid_region <- function(){
## Load the public linelist data by covid region ------------
##  (11 regions)
raw_dat_cr <- read_csv('../data/idph_public_covid_region.csv')  %>%
  rename(region = covid_region) %>%
  filter(region != 'unknown') %>%
  mutate(region = as.character(region)) %>%
  group_by(date, region) %>%
  summarise(cases = sum(cases)) %>%
  mutate(cases = ifelse(is.na(cases), 0, cases))  %>%
  ungroup
overall_dat_cr <- raw_dat_cr %>%
  group_by(date) %>%
  summarise(region = 'IL_Overall',
            cases = sum(cases))
bind_rows(raw_dat_cr, overall_dat_cr)%>%
  ungroup() %>%
  group_by(region) %>%
  arrange(date) %>%
  mutate(smoothed = smooth.spline(cases, spar = .6)$y %>% min_0,
         avg_7d = zoo::rollmean(cases, k = 7, fill = c(mean(cases[1:7], na.rm = T), 
                                                           NA, 
                                                           mean(cases[length(cases)-(0:6)], na.rm = T)))
         )%>%
  filter(date <= max(date)) %>%
  ungroup() %>%
  mutate(region = factor(region, levels = c(as.character(1:11), 'IL_Overall'))) %>%
  filter(date >= lubridate::as_date('2020-04-01'))%>%
  mutate(cases = ifelse(cases<0, 0, cases))
}

load_idph_public_adjusted_cases_covid_region <- function(){
  ## Load the public linelist data by covid region ------------
  ##  (11 regions)
  raw_dat_cr <- read_csv('../data/cbts_reg2_adjcases_noroll.csv')  %>%
    rename(region = covid_region) %>%
    filter(region != 'unknown') %>%
    mutate(region = as.character(region)) %>%
    group_by(date, region) %>%
    summarise(cases = sum(adjusted_cases)) %>%
    mutate(cases = ifelse(is.na(cases), 0, cases))  %>%
    ungroup
  overall_dat_cr <- raw_dat_cr %>%
    group_by(date) %>%
    summarise(region = 'IL_Overall',
              cases = sum(cases))
  bind_rows(raw_dat_cr, overall_dat_cr)%>%
    ungroup() %>%
    group_by(region) %>%
    arrange(date) %>%
    mutate(smoothed = smooth.spline(cases, spar = .6)$y %>% min_0,
           avg_7d = zoo::rollmean(cases, k = 7, fill = c(mean(cases[1:7], na.rm = T), 
                                                         NA, 
                                                         mean(cases[length(cases)-(0:6)], na.rm = T)))
    )%>%
    filter(date <= max(date)) %>%
    ungroup() %>%
    mutate(region = factor(region, levels = c(as.character(1:11), 'IL_Overall'))) %>%
    filter(date >= lubridate::as_date('2020-04-01'))%>%
    mutate(cases = ifelse(cases<0, 0, cases))
}


load_EPIC_admissions <- function(){
  ## Load the public linelist data by restore region ------------
  ## More up to date, public linelist
  read_csv('../data/cli_region11.csv') %>%
    mutate(date = format(as.Date(date), "%m/%d/%y")) %>%
    rename(
           nadmit = cli) %>%
    group_by(date) %>%
    summarise(nadmit = sum(nadmit)) %>%
    ungroup %>%
    tidyr::extract(date, into = c('month', 'day', 'year'), '(\\d\\d?)/(\\d\\d?)/(\\d\\d)') %>%
    mutate(nadmit = ifelse(is.na(nadmit), 0, nadmit),
           region = 'ALL_EPIC_HOSPITALS',
           date = as.Date(sprintf('%2i20-%2i-%2i', as.numeric(year), as.numeric(month), as.numeric(day)))) %>%
    mutate(smoothed = smooth.spline(nadmit, spar = .5)$y %>% min_0,
           avg_7d = zoo::rollmean(nadmit, k = 7, fill = c(mean(nadmit[1:7], na.rm = T), 
                                                             NA, 
                                                             mean(nadmit[length(nadmit)-(0:6)], na.rm = T)))
    )%>%
    select(-month,-day,-year) %>%
    ungroup() 
}


load_cli <- function(){
  ## Load the public linelist data by restore region ------------
  ## More up to date, public linelist
  region_cli = read.csv('../data/cli_admissions_latest.csv') %>%
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
