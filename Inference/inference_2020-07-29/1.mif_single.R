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
n_mif = 50
n_particles_mif = 3000
cooling_rate = 0.95
num_points = 100

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
region_to_test=as.numeric(args[5])
output_dir=args[6]


# Read in inference functions
source(covid_get_path(inference_file))
source('set_up_covariates_and_data.R')
pars$region_to_test = region_to_test

parnames = c('beta1_logit_', 'beta2_', 'num_init_', 'scale_phase3_')
lower_pars = c(0.5, 0.02, 10, 0.1)
upper_pars = c(0.9, 0.04, 1000, 0.5)
names(lower_pars) = paste0(parnames, region_to_test)
names(upper_pars) = paste0(parnames, region_to_test)
design=sobolDesign(lower=lower_pars, upper=upper_pars, num_points)

if (region_to_test == 1){
  rw_vec = rw.sd(
                 beta1_logit_1 = 0.01,
                 beta2_1 = 0.01,
                 scale_phase3_1=0.01,
                 num_init_1=ivp(0.1)
                 )

} else if(region_to_test == 2){
  rw_vec = rw.sd(
                 beta1_logit_2 = 0.01,
                 beta2_2 = 0.01,
                 scale_phase3_2=0.01,
                 num_init_2=ivp(0.1)
                 )

} else if (region_to_test == 3){
  rw_vec = rw.sd(
                 beta1_logit_3 = 0.01,
                 beta2_3 = 0.01,
                 scale_phase3_3=0.01,
                 num_init_3=ivp(0.1)
                 )
} else if (region_to_test == 4){
  rw_vec = rw.sd(
                 beta1_logit_4 = 0.01,
                 beta2_4 = 0.01,
                 scale_phase3_4=0.01,
                 num_init_4=ivp(0.1)
                 )

} else if (region_to_test == 5){
  rw_vec = rw.sd(
                 beta1_logit_5 = 0.01,
                 beta2_5 = 0.01,
                 scale_phase3_5=0.01,
                 num_init_5=ivp(0.1)
                 )
}

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

      print('pomp object created, running mifs')
      ## Run inference and return dataframes for output
      run_mif_and_output_result(
        jobid,
        arrayid,
        output_dir,
        pomp_object,
        dmeasFile,
        simstart,
        min_data_time,
        intervention_start,
        n_mif,
        rw_vec,
        cooling_rate,
        n_particles_mif)

      output_filename = sprintf("%s/output_mif_%s_%s.csv.rds", output_dir, jobid, arrayid)
      mf = readRDS(output_filename)


    ll <- logmeanexp(replicate(n_reps_pfilter,logLik(pfilter(mf,Np=n_particles_pfilter))),se=TRUE)
    finaloutput = (data.frame(as.list(mf@params),
                        loglik=ll[1],
                        loglik_se=ll[2],
                        model=dmeasFile,
                        nparticles=n_particles_pfilter,
                        n_reps=n_reps_pfilter))
    write.table(finaloutput, file=paste0(output_filename,'.csv'), sep = ",",col.names = TRUE, row.names=FALSE) 
    }
} ->dummy