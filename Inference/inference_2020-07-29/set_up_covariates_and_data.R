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
dmeasFile = covid_get_path(dmeasFile)
rmeasFile = covid_get_path(rmeasFile)

fraction_underreported_file = covid_get_path(fraction_underreported_file)

contact_filename = covid_get_path(contact_filename)

print('Reading in CSnippets')
# CSnippets

dmeasure_snippet <- read_Csnippet(dmeasFile)
rprocess_snippet <- read_Csnippet(rprocFile)
rinit_snippet <- read_Csnippet(initFile)
rmeasure_snippet <- read_Csnippet(rmeasFile)

print('Setting beta scaling')
beta_scales = get_scale(t_logistic_start=simend+100,
                        intervention_lift=simend+100,
                        simstart=simstart,
                        simend=simend,
                        max_scales=c(1,1,1,1,1))
    
print('Getting fraction underreported')
fraction_underreported = read.csv(fraction_underreported_file)
row.names(fraction_underreported) = fraction_underreported$time
# icu reporting


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
population5 = read.csv(covid_get_path(population_filename_5))
colnames(population5) <- c("AGE_GROUP_MIN", "POPULATION")

n_age_groups = nrow(population1)

print('Initializing parameters')
### User-specification of parameters and interventions for model 
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)

age_dist_frame = read.csv(covid_get_path(age_dist_file)) 

print('Adding intervention scaling')
pars = add_interventions(covid_get_path(intervention_file), pars)

print('Loading data')
# Load real data, assume read in from civis

if (pars$n_regions == 4){
        region_order = c('northcentral','central','northeast','southern')
        population_list=list(population1, population2, population3, population4)
    } else if (pars$n_regions == 5){
       region_order = c('northcentral','central','northeast','southern', 'chicago')
       population_list=list(population1, population2, population3, population4, population5)
    }

civis_data = read.csv(data_filename)  %>% 
    mutate(Date=as.Date(date), total_deaths = hosp_deaths+nonhosp_deaths)

get_idph = function(date, region){
    idph[which(idph$Date==date & idph$restore_region==region), 'total_deaths']
}

# Add idph total deaths to total deaths
idph = read.csv(idph_filename) %>% 
    mutate(Date=as.Date(test_date), restore_region=region, total_deaths=inc_deaths) %>%
    select(Date, restore_region, total_deaths)

df = civis_data %>% mutate(
    source=case_when((is.na(total_deaths)) ~ 'Public',
                     (!is.na(total_deaths))~ 'IDPH line list'),
    total_deaths=case_when((is.na(total_deaths)) ~ as.numeric(mapply(get_idph, Date, restore_region)),
                           (!is.na(total_deaths))~ as.numeric(total_deaths))) %>%
    mutate(time = as.numeric(Date - as.Date(t_ref))) %>%
    group_by(restore_region) %>%
    mutate(confirmed_covid_icu = predict(loess(confirmed_covid_icu ~ time, span=0.75), time),
           total_deaths = predict(loess(total_deaths ~ time, span=0.75), time),
           covid_non_icu = predict(loess(covid_non_icu ~ time, span=0.75), time),
           date=Date) %>%
    ungroup() %>%
    select(restore_region, confirmed_covid_icu, covid_non_icu, total_deaths, date) %>%
    mutate(confirmed_covid_icu = if_else(confirmed_covid_icu < 0, 0, round(confirmed_covid_icu)),
           total_deaths = if_else(total_deaths < 0, 0, round(total_deaths)),
           covid_non_icu = if_else(covid_non_icu < 0, 0, round(covid_non_icu)))


df_ICU = df %>% select(date, restore_region, confirmed_covid_icu) %>% 
    spread(restore_region, confirmed_covid_icu) %>% 
    select(date, region_order) %>% 
    filter(date>=min_data_time_ICU)
names(df_ICU) = c('time', paste0('ObsICU_', 1:pars$n_regions))

df_death = df %>% 
    select(date, restore_region, total_deaths) %>% 
    spread(restore_region, total_deaths) %>% 
    select(date, region_order)
names(df_death) = c('time', paste0('ObsDeaths_', 1:pars$n_regions))

df_hosp = df %>% 
    select(date, restore_region, covid_non_icu) %>% 
    spread(restore_region, covid_non_icu) %>% 
    select(date, region_order)
names(df_hosp) = c('time', paste0('ObsHosp_', 1:pars$n_regions))

data = left_join(df_death, df_ICU, by='time')
data = left_join(data, df_hosp, by='time')
data$time = as.numeric(as.Date(data$time) - as.Date(t_ref))
print(tail(data))

print('Set up parameter transformation')
observed_names=names(data %>% select(-time))
transformation=parameter_trans(
    log=c('num_init_1',
    'num_init_2',
    'num_init_3',
    'num_init_4',
    'num_init_5',
    'beta2_1',
    'beta2_2',
    'beta2_3',
    'beta2_4',
    'beta2_5',
    'scale_phase3_1',
    'scale_phase3_2',
    'scale_phase3_3',
    'scale_phase3_4',
    'scale_phase3_5'),

    logit=c('beta1_logit_1',
    'beta1_logit_2',
    'beta1_logit_3',
    'beta1_logit_4',
    'beta1_logit_5')
    )


