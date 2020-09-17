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

cooling_rate = 0.99
num_points = 100
n_reps_pfilter = 3


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
print('Finished covariate setup')


if (use_changepoint){
  parnames = c('beta1_', 'beta2_', 'num_init_', 'scale_phase3_', 'scale_phase4_', 'scale_phase5_')
  final_parnames = paste0(parnames, region_to_test)
  design = read.csv('final_mle_pars.csv') %>% filter(param_name %in% final_parnames)
  temp = unlist(design$value, use.names=F)
  names(temp) = unlist(design$param_name, use.names=F)
  temp = t(temp)
  design = data.frame(temp) %>% 
    slice(rep(1, each = num_points)) 

  #design = readRDS('2020-08-26.final.mle_grid.rds') %>% 
  #  filter(region_to_test==region_to_test) %>% 
  #  slice(rep(1, each = num_points)) %>% 
  #  select(-region_to_test)
  #names(design) = paste0(parnames, region_to_test)
  print(design)

  if (region_to_test == 1){
    rw_vec = rw.sd(
                   beta1_1 = ifelse(time < 62, 0.007, 0),
                   beta2_1 = ifelse(time >= 62 & time < 85, 0.007, 0),
                   scale_phase3_1=ifelse(time >= 85 & time < 151, 0.007, 0),
                   scale_phase4_1=ifelse(time >= 151 & time < 200, 0.007, 0),
                   scale_phase5_1=ifelse(time >= 200, 0.007, 0),                
                   num_init_1=ivp(0.05)
                   )

  } else if(region_to_test == 2){
    rw_vec = rw.sd(
                   beta1_2 = ifelse(time < 62, 0.007, 0),
                   beta2_2 = ifelse(time >= 62 & time < 85, 0.007, 0),
                   scale_phase3_2=ifelse(time >= 85 & time < 151, 0.007, 0),
                   scale_phase4_2=ifelse(time >= 151 & time < 200, 0.007, 0),
                   scale_phase5_2=ifelse(time >= 200, 0.007, 0),
                   num_init_2=ivp(0.05)
                   )

  } else if (region_to_test == 3){
    rw_vec = rw.sd(
                   beta1_3 = ifelse(time < 62, 0.007, 0),
                   beta2_3 = ifelse(time >= 62 & time < 85, 0.007, 0),
                   scale_phase3_3=ifelse(time >= 85 & time < 151, 0.007, 0),
                   scale_phase4_3=ifelse(time >= 151 & time < 200, 0.007, 0),
                   scale_phase5_3=ifelse(time >= 200, 0.007, 0),   
                   num_init_3=ivp(0.05)
                   )
  } else if (region_to_test == 4){
    rw_vec = rw.sd(
                   beta1_4 = ifelse(time < 62, 0.007, 0),
                   beta2_4 = ifelse(time >= 62 & time < 92, 0.007, 0),
                   scale_phase3_4=ifelse(time >= 92 & time < 151, 0.007, 0),
                   scale_phase4_4=ifelse(time >= 151 & time < 200, 0.007, 0),
                   scale_phase5_4=ifelse(time >= 200, 0.007, 0),  
                   num_init_4=ivp(0.05)
                   )

  } else if (region_to_test == 5){
    rw_vec = rw.sd(
                   beta1_5 = ifelse(time < 62, 0.007, 0),
                   beta2_5 = ifelse(time >= 62 & time < 85, 0.007, 0),
                   scale_phase3_5=ifelse(time >= 85 & time < 151, 0.007, 0),
                   scale_phase4_5=ifelse(time >= 151 & time < 200, 0.007, 0),
                   scale_phase5_5=ifelse(time >= 200, 0.007, 0),   
                   num_init_5=ivp(0.05)
                   )
  }else if (region_to_test == 6){
    rw_vec = rw.sd(
                   beta1_6 = ifelse(time < 62, 0.007, 0),
                   beta2_6 = ifelse(time >= 62 & time < 92, 0.007, 0),
                   scale_phase3_6=ifelse(time >= 92 & time < 151, 0.007, 0),
                   scale_phase4_6=ifelse(time >= 151 & time < 200, 0.007, 0),
                   scale_phase5_6=ifelse(time >= 200, 0.007, 0),   
                   num_init_6=ivp(0.05)
                   )
}
}

cl <- makeCluster(num_cores)
registerDoParallel(cl)

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
        pars[[name]] = abs(rnorm(1, design[jobid, name], design[jobid, name] * 0.02))
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
          beta_covar=beta_covariate_column,
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