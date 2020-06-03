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

args = commandArgs(trailingOnly=TRUE)

root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

default_par_file = './default_parameter_values.csv'
deltaT = 0.1

jobid_master = as.numeric(args[1])
arrayid = as.numeric(args[2])
output_dir = args[3]
dmeasFile = args[4]
data_filename = args[5]
population_filename_1 = args[6]
population_filename_2 = args[7]
population_filename_3 = args[8]
simstart = args[9]
simend = args[10]
min_data_time = args[11]
intervention_start=args[12]
num_cores=as.numeric(args[13])
maxjobs=as.numeric(args[14])
min_data_time_ICU = args[15]
init_file=args[16]
rprocFile=args[17]
nu_scales_file=args[18]
fraction_underreported_file=args[19]
contact_filename=args[20]
inference_file=args[21]
t_ref=args[22]
intervention_file=args[23]
population_filename_4 = args[24]
region_to_test=as.numeric(args[25])

n_reps_pfilter=5
n_particles_pfilter=3000

design = read.csv('points_to_search.csv')

num_points=nrow(design)
print(head(design))
source(covid_get_path(inference_file))

cl <- makeCluster(num_cores)
registerDoParallel(cl)

source('set_up_covariates_and_data.R')

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

    if (jobid <= num_points){
      for (name in names(design)){
        pars[[name]] = design[jobid, name]
      }
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
          obsnames=observed_names,
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