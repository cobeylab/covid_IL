## Run mif searches, can be implemented in parallel on the cluster 
## Sylvia Ranjeva, Phil Arevalo
## April 2020
## -----------------------------------------------------

## Load packages and specify files 
library(dplyr)
library(pomp)
library(ggplot2)
library(magrittr)
library(lubridate)
library(MASS)
library(reshape2)
library(tidyr)
library(foreach)
library(doParallel)

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

# Set mif search parameters
n_reps_pfilter = 3
n_particles_pfilter = 5000

root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

source('./input_file_specification.R')
default_par_file = './default_parameter_values.csv'

# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1

# Reading input args
args = commandArgs(trailingOnly=TRUE)
jobid_master = as.numeric(args[1])
arrayid = as.numeric(args[2])
num_cores=as.numeric(args[3])
maxjobs=as.numeric(args[4])
output_dir=args[5]


# Read in inference functions
source(covid_get_path(inference_file))
source('set_up_covariates_and_data.R')

design=readRDS('mle_grid.rds')
num_points = nrow(design)


cl <- makeCluster(num_cores)
registerDoParallel(cl)


# Setting up intervention
pars = add_interventions(covid_get_path(intervention_file), pars)

temp_age_dist = age_dist_frame %>% 
  filter(b_elderly == pars$b_elderly) %>% select(-b_elderly)
age_dist = as.numeric(temp_age_dist$value)
names(age_dist) = temp_age_dist$param_name
pars = c(age_dist, pars)

print(pars)

loopstart = (jobid_master - 1) * maxjobs + 1
loopend = loopstart + (maxjobs-1)
print(loopstart)
print(loopend)
print(maxjobs)
print(num_cores)
foreach(jobid=loopstart:loopend,
  .packages=c('pomp', 'dplyr', 'reshape2')) %dopar%{
    ## Set parameters to search over
    select <- dplyr::select
    rename <- dplyr::rename
    summarize <- dplyr::summarise
    contains <- dplyr::contains

    if (jobid <= num_points){
      for (name in names(design)){
        if(name == 'region_to_test'){
          pars[['region_to_test']] = design[jobid, name]
        } else{
          pars[[paste0(name, '_', design[jobid, 'region_to_test'])]] = design[jobid, name]
        }
        
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

      print('pomp object created, running grid')

    ll <- logmeanexp(replicate(n_reps_pfilter,logLik(pfilter(pomp_object,Np=n_particles_pfilter))),se=TRUE)
    finaloutput = (data.frame(as.list(pomp_object@params),
                        loglik=ll[1],
                        loglik_se=ll[2],
                        model=dmeasFile,
                        nparticles=n_particles_pfilter,
                        n_reps=n_reps_pfilter))
    write.table(finaloutput, file=sprintf("%s/output_grid_%s_%s.csv", output_dir, jobid, arrayid), sep = ",",col.names = TRUE, row.names=FALSE) 
    }
} ->dummy