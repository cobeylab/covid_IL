# Create result directory if it doesn't exist
dir.create(file.path(output_dir), showWarnings = FALSE)
human_output_dir = sprintf("./%s_outputs/", model_name)
dir.create(file.path(human_output_dir), showWarnings = FALSE)

## Setting dates
print('Setting dates')
simstart = convert_date_to_time(simstart, t_ref)
simend = convert_date_to_time(simend, t_ref)
min_data_time = convert_date_to_time(min_data_time, t_ref)
intervention_start = convert_date_to_time(intervention_start, t_ref)
project_end = convert_date_to_time(project_end, t_ref)
min_data_time_ICU = as.Date(min_data_time_ICU)

## Filenames
print('Reading in files')
initFile = covid_get_path(init_file)
rprocFile = covid_get_path(rprocFile)
dmeasFile = covid_get_path(dmeasFile)
rmeasFile = covid_get_path(rmeasFile)

# covariates
fraction_underreported_filename = covid_get_path(fraction_underreported_file)
nonhosp_death_filename = covid_get_path(nonhosp_deaths)
beta_covariate_filename = covid_get_path(beta_covariate)
beta_scale_filename = covid_get_path(beta_scale)
emr_reporting_filename = covid_get_path(emr_report)

## CSnippets
print('Reading in CSnippets')
dmeasure_snippet <- read_Csnippet(dmeasFile)
rprocess_snippet <- read_Csnippet(rprocFile)
rinit_snippet <- read_Csnippet(initFile)
rmeasure_snippet <- read_Csnippet(rmeasFile)
    
## Fraction underreported
print('Loading covariates')
fraction_underreported = read.csv(fraction_underreported_filename)
nonhosp_deaths = read.csv(nonhosp_death_filename)
beta_cov = read.csv(beta_covariate_filename)
beta_scale_covar = read.csv(beta_scale_filename)
emr_report_covar = read.csv(emr_reporting_filename) %>%
    filter(covid_region %in% as.character(1:11)) %>%
    select(time, covid_region, death_reporting, hosp_reporting) %>%
    gather(Compartment, Reporting, death_reporting, hosp_reporting) %>%
    mutate(Compartment = paste0(Compartment, '_', covid_region)) %>%
    select(-covid_region) %>%
    spread(Compartment, Reporting) %>%
    select(time, paste0('hosp_reporting_', 1:11), paste0('death_reporting_', 1:11))
    
## Population sizes
print('Reading in population info')
population_list=read.csv(covid_get_path(population_file))
n_age_groups = 1
stopifnot(n_age_groups == 1)

## User-specification of parameters
print('Initializing parameters')
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)
pars$tmax = simend
pars$tmin = simstart

## read in fitting data
hosp = read.csv(covid_get_path(emr_data_file)) %>%
  filter(covid_region %in% 1:11) %>%
  select(time, covid_region, ObsHosp_1) %>%
  spread(covid_region, ObsHosp_1) %>%
  select(time, as.character(1:11))
names(hosp) = c('time', paste0('ObsHosp_', 1:11)) 
hosp_deaths = read.csv(covid_get_path(emr_data_file)) %>%
  filter(covid_region %in% 1:11) %>%
  select(time, covid_region, ObsHospDeaths_1) %>%
  spread(covid_region, ObsHospDeaths_1) %>%
  select(time, as.character(1:11))
names(hosp_deaths) = c('time', paste0('ObsHospDeaths_', 1:11))
ll_deaths = read.csv(covid_get_path(ll_file)) %>%
  mutate(time = as.numeric(as.Date(date) - as.Date('2020-01-14'))) %>%
  select(time, covid_region, ll_deaths) %>%
  spread(covid_region, ll_deaths) %>%
  select(time, as.character(1:11))
names(ll_deaths) = c('time', paste0('ObsDeaths_', 1:11))
fitting_data = full_join(hosp, hosp_deaths, by='time') %>%
  full_join(ll_deaths, by='time') %>%
  arrange(time)
NA_row = data.frame(time=min(fitting_data$time) -1)
fitting_data = full_join(NA_row, fitting_data)