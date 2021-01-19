
## Run mif searches, can be implemented in parallel on the cluster 
## Sylvia Ranjeva, Phil Arevalo
## April 2020
## -----------------------------------------------------

## Load packages and specify files 
library(dplyr)
library(pomp)
library(tidyr)
library(foreach)
library(doParallel)

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)
# Read in inference functions
source('./input_file_specification_no_age.R')
source(covid_get_path(function_file))
source('./set_up_covariates_and_data_no_age.R')


# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))

deltaT = 0.1

# Reading input args
args = commandArgs(trailingOnly=TRUE)
arrayid = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
num_cores = as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
region_to_test = arrayid

stopifnot(region_to_test %in% c(-1, 1:11))

pars$region_to_test = region_to_test
regtemp = region_to_test
pars$tmin = simstart
sim_times = data.frame(time=47:project_end)
pars$tmax = project_end

sim_times = data.frame(time=47:project_end)
fitting_data = fitting_data %>%
  full_join(sim_times, by='time') %>%
  arrange(time)
finalobsnames = names(fitting_data %>% select(-time))
max_report_time = max(emr_report_covar$time)
covar_region =  emr_report_covar %>%
  select(time, all_of(grep(pattern=sprintf("_%s$", region_to_test), names(.)))) %>%
  full_join(sim_times, by='time') %>%
  arrange(time)
covar_future = covar_region %>%
  filter(time > max_report_time) %>%
  mutate_at(vars(-time), ~0.95)
covar_region = covar_region %>%
  filter(time <= max_report_time) %>%
  bind_rows(covar_future) %>%
  arrange(time)

files = Sys.glob(sprintf('%s/sim_pfilter_reg_%s_*.csv', output_dir, region_to_test))

registerDoParallel(cores=num_cores)
foreach(i=1:length(files), .combine='rbind') %dopar%{
  pf.new = readRDS(files[i])
  traj = filter.traj(pf.new)
  times = as.numeric(names(traj[1,1,])) # get the times
  data.frame(t(traj[,1,])) %>%
      select(all_of(grep(pattern=sprintf("_%s$", region_to_test), names(.)))) %>%
      mutate(time = times,
             .id = i)
} -> trajectory

print('Outputting results')

hosp_census_cols =  names(trajectory)[grep(pattern="^IH[1,2]_.*", names(trajectory))]
icu_census_cols = names(trajectory)[grep(pattern="^IC_.*", names(trajectory))]
new_death_cols = names(trajectory)[grep(pattern="^new_deaths_.*", names(trajectory))]
new_hosp_death_cols = names(trajectory)[grep(pattern="^new_hosp_deaths_.*", names(trajectory))]

obscolnames = c('ObsHosp','ObsICU','ObsHospDeaths','ObsDeaths')
newnames = setNames(obscolnames, paste0(obscolnames, '_', region_to_test))
sim_full = trajectory %>%
  left_join(covar_region, by = 'time') %>% # join in covariates
  mutate(HospCensus = trajectory %>% select(all_of(hosp_census_cols)) %>% rowSums(),
    ICUCensus = trajectory %>% select(all_of(icu_census_cols)) %>% rowSums(),
    Deaths = trajectory %>% select(all_of(new_death_cols)) %>% rowSums(),
    HospDeaths = trajectory %>% select(all_of(new_hosp_death_cols)) %>% rowSums()) %>% # get the total counts across subcompartments etc.
  mutate(ObsHosp = mapply(FUN=rpois, 1, !!as.name(paste0('hosp_reporting_', region_to_test)) * (HospCensus)),
    ObsICU = mapply(FUN=rpois, 1, !!as.name(paste0('icu_reporting_', region_to_test)) * (ICUCensus)),
    ObsHospDeaths = mapply(FUN=rpois, 1, 0.99 * (HospDeaths)),
    ObsDeaths = mapply(FUN=rpois, 1, !!as.name(paste0('DeathReportTrack_1_', region_to_test)) * (Deaths))) %>%
  rename_(.dots = newnames) %>%
  select(-one_of(c(covar_region %>% select(-time) %>% names(), 'HospCensus', 'ICUCensus', 'Deaths', 'HospDeaths')))

write.csv(sim_full, sprintf("%s/full_projection_%s_%s.csv", human_output_dir, model_name, region_to_test), row.names=F)
