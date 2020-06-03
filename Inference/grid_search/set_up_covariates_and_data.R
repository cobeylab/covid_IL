library(tidyverse)

# Create result directory if it doesn't exist
dir.create(file.path(output_dir), showWarnings = FALSE)

# Setting dates
print('Setting dates')
simstart = convert_date_to_time(simstart, t_ref)
simend = convert_date_to_time(simend, t_ref)
min_data_time = convert_date_to_time(min_data_time, t_ref)
intervention_start = convert_date_to_time(intervention_start, t_ref)
min_data_time_ICU = as.Date(min_data_time_ICU)


print('Reading in files')
initFile = covid_get_path(init_file)
rprocFile = covid_get_path(rprocFile)
nu_scales_file = covid_get_path(nu_scales_file)
fraction_underreported_file = covid_get_path(fraction_underreported_file)
contact_filename = covid_get_path(contact_filename)

print('Reading in CSnippets')
# CSnippets
dmeasFile = covid_get_path(dmeasFile)
dmeasure_snippet <- read_Csnippet(dmeasFile)
rprocess_snippet <- read_Csnippet(rprocFile)
rinit_snippet <- read_Csnippet(initFile)
rmeasure_snippet <- read_Csnippet(covid_get_path('Inference/rmeasure_hosp_non_hosp_ICU.c'))

print('Reading in nu scaling')
# import scaling for nu
nu_scales = read.csv(nu_scales_file)
row.names(nu_scales) = nu_scales$time

print('Setting beta scaling')
beta_scales = get_scale(t_logistic_start=simend+100,
                        intervention_lift=simend+100,
                        simstart=simstart,
                        simend=simend,
                        max_scales=c(1,1,1,1))
    
print('Getting fraction underreported')
fraction_underreported = read.csv(fraction_underreported_file)
row.names(fraction_underreported) = fraction_underreported$time

print('Loading contact matrix')
load(contact_filename)


print('Reading in population info')
population1 = read.csv(covid_get_path(population_filename_1))
colnames(population1) <- c("AGE_GROUP_MIN", "POPULATION")
population2 = read.csv(covid_get_path(population_filename_2))
colnames(population2) <- c("AGE_GROUP_MIN", "POPULATION")
population3 = read.csv(covid_get_path(population_filename_3))
colnames(population3) <- c("AGE_GROUP_MIN", "POPULATION")
population4 = read.csv(covid_get_path(population_filename_4))
colnames(population4) <- c("AGE_GROUP_MIN", "POPULATION")
population_list=list(population1, population2, population3, population4)

n_age_groups = nrow(population1)

print('Loading data')
# Load real data, assume read in from civis
region_order = c('northcentral','central','northeast','southern')
df = read.csv(data_filename) %>% mutate(date = as.Date(date))
df_ICU = df %>% select(date, restore_region, confirmed_covid_icu) %>% spread(restore_region, confirmed_covid_icu) %>% select(date, region_order) %>% filter(date>=min_data_time_ICU)
names(df_ICU) = c('time', paste0('ObsICU_', 1:4))
df_hosp = df %>% select(date, restore_region, hosp_deaths) %>% spread(restore_region, hosp_deaths) %>% select(date, region_order)
names(df_hosp) = c('time', paste0('ObsHospDeaths_', 1:4))
df_nonhosp = df %>% select(date, restore_region, nonhosp_deaths) %>% spread(restore_region, nonhosp_deaths) %>% select(date, region_order)
names(df_nonhosp) = c('time', paste0('ObsNonHospDeaths_', 1:4))
data = left_join(df_hosp, df_ICU, 'time') %>% inner_join(df_nonhosp, 'time')
data$time = as.numeric(as.Date(data$time) - as.Date(t_ref))
print(data)

print('Initializing parameters')
### User-specification of parameters and interventions for model 
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)
pars$region_to_test = region_to_test


print('Adding intervention scaling')
pars = add_interventions(covid_get_path(intervention_file), pars)

print('Set up parameter transformation')
observed_names=names(data %>% select(-time))
transformation=parameter_trans(log=c('beta1_1',
    'beta1_2',
    'beta1_3',
    'beta1_4',
    'beta2_1',
    'beta2_2',
    'beta2_3',
    'beta2_4',
    'num_init_1',
    'num_init_2',
    'num_init_3',
    'num_init_4'),
    logit=c('region_non_hosp_1',
        'region_non_hosp_2',
        'region_non_hosp_3',
        'region_non_hosp_4'))


