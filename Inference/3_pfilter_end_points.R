## Run mif searches, can be implemented in parallel on the cluster 
## Sylvia Ranjeva
## April 2020
## -----------------------------------------------------

## Load packages and specify files 
require(dplyr)
require(pomp)
require(ggplot2)
require(magrittr)
require(lubridate)
require(MASS)
require(reshape2)
require(dplyr)
require(foreach)
require(doParallel)


select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

# pfilter parameters
n_reps_pfilter=5
n_particles_pfilter=10000

stopifnot(file.exists('mif_50iter_end_points.csv'))
design = read.csv('mif_50iter_end_points.csv')

root <- '../../../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)
default_par_file = './default_parameter_values.csv'
# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))
deltaT = 0.1

# Read in input args
args = commandArgs(trailingOnly=TRUE)
jobid_master = as.numeric(args[1])
arrayid = as.numeric(args[2])
output_dir = args[3]
dmeasFile = args[4]
data_filename_1 = args[5]
data_filename_2 = args[6]
data_filename_3 = args[7]
population_filename_1 = args[8]
population_filename_2 = args[9]
population_filename_3 = args[10]
simstart = args[11]
simend = args[12]
min_data_time = args[13]
intervention_start=args[14]
num_cores=as.numeric(args[15])
maxjobs=as.numeric(args[16])
min_data_time_ICU = args[17]
init_file=args[18]
rprocFile=args[19]
nu_scales_file=args[20]
fraction_underreported_file=args[21]
contact_filename=args[22]
inference_file=args[23]
t_ref=args[24]
intervention_file=args[25]
source(covid_get_path(inference_file))

# Setting dates
simstart = convert_date_to_time(simstart, t_ref)
simend = convert_date_to_time(simend, t_ref)
min_data_time = convert_date_to_time(min_data_time, t_ref)
intervention_start = convert_date_to_time(intervention_start, t_ref)
min_data_time_ICU = convert_date_to_time(min_data_time_ICU, t_ref)

cl <- makeCluster(num_cores)
registerDoParallel(cl)

print('Reading in files')
initFile = covid_get_path(init_file)
rprocFile = covid_get_path(rprocFile)
nu_scales_file = covid_get_path(nu_scales_file)
fraction_underreported_file = covid_get_path(fraction_underreported_file)
contact_filename = covid_get_path(contact_filename)

# CSnippets
dmeasFile = covid_get_path(dmeasFile)
dmeasure_snippet <- read_Csnippet(dmeasFile)
rprocess_snippet <- read_Csnippet(rprocFile)
rinit_snippet <- read_Csnippet(initFile)

# import scaling for nu
nu_scales = read.csv(nu_scales_file)
row.names(nu_scales) = nu_scales$time

# Set beta scales to be the same as nu scales, i.e., 1 all over
beta_scales = get_scale(simend+100,
                        simend+100,
                        simstart,
                        simend,
                        1)
    
# import fraction underreported covariate
fraction_underreported = read.csv(fraction_underreported_file)
row.names(fraction_underreported) = fraction_underreported$time

# Import contact matrix data
load(contact_filename)

# Make population a list that has the regions concatenated together
population1 = read.csv(covid_get_path(population_filename_1))
colnames(population1) <- c("AGE_GROUP_MIN", "POPULATION")
population2 = read.csv(covid_get_path(population_filename_2))
colnames(population2) <- c("AGE_GROUP_MIN", "POPULATION")
population3 = read.csv(covid_get_path(population_filename_3))
colnames(population3) <- c("AGE_GROUP_MIN", "POPULATION")
population_list=list(population1, population2, population3)
n_age_groups = nrow(population1)

# Load real data, assume read in from civis
data_1 = get_civis_data(data_filename_1)  %>% mutate(ObsDeaths_1=ObsDeaths, ObsICU_1=ObsICU) %>% select(time, ObsDeaths_1, ObsICU_1)
data_2 = get_civis_data(data_filename_2)  %>% mutate(ObsDeaths_2=ObsDeaths, ObsICU_2=ObsICU) %>% select(time, ObsDeaths_2, ObsICU_2)
data_3 = get_civis_data(data_filename_3)  %>% mutate(ObsDeaths_3=ObsDeaths, ObsICU_3=ObsICU) %>% select(time, ObsDeaths_3, ObsICU_3)
data = merge(data_1, data_2, by='time')
data = merge(data, data_3, by='time')
# Add NAs if necessary
data <- data %>% filter(time>=simstart, time<=simend) %>% 
                 mutate(
                        ObsDeaths_1=ifelse(time<min_data_time, NA, ObsDeaths_1),
                        ObsDeaths_2=ifelse(time<min_data_time, NA, ObsDeaths_2),
                        ObsDeaths_3=ifelse(time<min_data_time, NA, ObsDeaths_3),

                        ObsICU_1=ifelse(time<min_data_time_ICU, NA, ObsICU_1),
                        ObsICU_2=ifelse(time<min_data_time_ICU, NA, ObsICU_2),
                        ObsICU_3=ifelse(time<min_data_time_ICU, NA, ObsICU_3)) %>%
                select(ObsDeaths_1, ObsDeaths_2, ObsDeaths_3, ObsICU_1, ObsICU_2, ObsICU_3, time)
print(data)

print('Initializing parameters')
### User-specification of parameters and interventions for model 
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)

# Add in intervention scaling
pars = add_interventions(covid_get_path(intervention_file), pars)
print(pars)

loopstart = (jobid_master - 1) * maxjobs + 1
loopend = loopstart + (maxjobs-1)
print(loopstart)
print(loopend)

foreach(jobid=loopstart:loopend,
  .packages=c('pomp', 'dplyr', 'reshape2')) %dopar%{
    ## Set parameters to search over
    select <- dplyr::select
    rename <- dplyr::rename
    summarize <- dplyr::summarise
    contains <- dplyr::contains

    if (jobid <= nrow(design)){
      pars$beta1 = design[jobid, 'beta1']
      pars$beta2_1 = design[jobid, 'beta2_1']
      pars$beta2_2 = design[jobid, 'beta2_2']
      pars$beta2_3 = design[jobid, 'beta2_3']
      pars$num_init_1 = round(design[jobid, 'num_init_1'])
      pars$num_init_2 = round(design[jobid, 'num_init_2'])
      pars$num_init_3 = round(design[jobid, 'num_init_3'])
      transformation=parameter_trans(log=c('beta1', 'beta2_1', 'beta2_2', 'beta2_3', 'num_init_1', 'num_init_2', 'num_init_3'))
      ## Make a pomp object for inference
      print('Making pomp object')
      pomp_object <- make_pomp_object_covid(
          n_regions = pars$n_regions,
          n_age_groups = n_age_groups,
          input_params = pars,
          delta_t = deltaT,
          contacts=pomp_contacts,
          population=population_list,
          nu_scales=nu_scales,
          beta_scales=beta_scales,
          frac_underreported=fraction_underreported,
          dmeasure_Csnippet = dmeasure_snippet,
          rprocess_Csnippet = rprocess_snippet,
          rinit_Csnippet = rinit_snippet,
          data=data,
          fitstart=simstart,
          time_column='time',
          obsnames=c('ObsDeaths_1', 'ObsDeaths_2', 'ObsDeaths_3', 'ObsICU_1', 'ObsICU_2', 'ObsICU_3'),
          transformations=transformation,
          inherit_parameters=FALSE
          )

      print('pomp object created, running pfilter')

      ## Run inference and return dataframes for output
      run_pfilter_and_output_result(n_reps_pfilter, 
        n_particles_pfilter, 
        jobid,
        arrayid,
        output_dir,
        pomp_object,
        dmeasFile,
        simstart,
        min_data_time,
        intervention_start)   
    }
     
} -> dummy_var