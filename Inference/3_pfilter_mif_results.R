## Run pfilters on mif objects
## Sylvia Ranjeva, Phil Arevalo
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


select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

args = commandArgs(trailingOnly=TRUE)


n_reps_pfilter=10
n_particles_pfilter=6000

root <- '../../../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

source('./inference_functions.R')

jobid = as.numeric(args[1])
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


dmeasFile = covid_get_path(dmeasFile)
file_list <- list.files("./results/", pattern = "*.rds$")
miffile = file_list[jobid]
mf <- readRDS(paste0('./results/', miffile))

## Run inference and return dataframes for output
run_pfilter_and_output_result(n_reps_pfilter, 
  n_particles_pfilter, 
  jobid,
  arrayid,
  output_dir,
  mf,
  miffile,
  simstart,
  min_data_time,
  intervention_start)       
