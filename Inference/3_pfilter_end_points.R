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

design = read.csv('mif_50iter_best_points.csv')
args = commandArgs(trailingOnly=TRUE)

n_reps_pfilter=5
n_particles_pfilter=5000

root <- '../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

source(covid_get_path('Inference/inference_functions.R'))
initFile = covid_get_path('Inference/initializer_compartment_distribute_IH4.c')
rprocFile = covid_get_path('Inference/rprocess_interventionbeta_IH4.c')
nu_scales_file = covid_get_path('Data/nu_scaling.csv')
fraction_underreported_file = covid_get_path('Data/frac_underreported.csv')
default_par_file = covid_get_path('Parameters/parameter_values.csv')
contact_filename = covid_get_path('Data/formatted_contacts_IL.RData')
deltaT = 0.1


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

simstart = as.numeric(args[11])
simend = as.numeric(args[12])
min_data_time = as.numeric(args[13])
intervention_start=as.numeric(args[14])

scalework=as.numeric(args[15])
scaleschool=as.numeric(args[16])
scalehome=as.numeric(args[17])
scaleother=as.numeric(args[18])
num_cores=as.numeric(args[19])
maxjobs=as.numeric(args[20])

min_data_time_ICU = as.numeric(args[21])

cl <- makeCluster(num_cores)
registerDoParallel(cl)

print('Reading in files')
# CSnippets
dmeasFile = covid_get_path(dmeasFile)
dmeasure_snippet <- read_Csnippet(dmeasFile)
rprocess_snippet <- read_Csnippet(rprocFile)
rinit_snippet <- read_Csnippet(initFile)

# import scaling for nu
nu_scales = read.csv(nu_scales_file)
row.names(nu_scales) = nu_scales$time

# import fraction underreported covariate
fraction_underreported = read.csv(fraction_underreported_file)
row.names(fraction_underreported) = fraction_underreported$time

# Import contact matrix data
load(contact_filename)

# Make population a vector that has the regions concatenated together
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

## Set intervention start
pars$t_start = intervention_start

## Set intervention strength
pars$scalework =scalework
pars$scaleschool=scaleschool
pars$scaleother=scaleother
pars$scalehome=scalehome

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