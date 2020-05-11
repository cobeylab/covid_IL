## Run mif searches, can be implemented in parallel on the cluster 
## Sylvia Ranjeva
## April 2020
## -----------------------------------------------------

## Load packages and specify files 
library(dplyr)
library(pomp)
library(ggplot2)
library(MASS)
library(reshape2)
library(dplyr)
library(rmutil)
library(tidyverse)
library(data.table)
library(foreach)

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

root <- '../../../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

default_par_file = './default_parameter_values.csv'
# Require the default parameter file to be in the same directory as the script
stopifnot(file.exists(default_par_file))
deltaT = 0.1

# Read in arguments
args = commandArgs(trailingOnly=TRUE)
output_dir = args[1]
dmeasFile = args[2]
data_filename_1 = args[3]
data_filename_2 = args[4]
data_filename_3 = args[5]
population_filename_1 = args[6]
population_filename_2 = args[7]
population_filename_3 = args[8]
simstart = args[9]
simend = args[10]
min_data_time = args[11]
intervention_start= args[12]
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
source(covid_get_path(inference_file))

# Aggregate points
print('Aggregating points')
output_filename <- "consolidated.mif.pfilter.csv"
liksample_filename <- "consolidated.mif.sample_from_surface.csv"
#file_path = output_dir #Replace this with a file path to the directory storing the csv files 
#setwd(file_path)
#file_list <- list.files(pattern = "*.csv$")
#data <- rbindlist(lapply(file_list, fread))
#df_output <- as.data.frame(data) %>% arrange(desc(loglik))
#design = df_output
#write.csv(df_output, file = paste0('../', output_filename))
#df_output$lik = exp(df_output$loglik - max(df_output$loglik))
#num_sims = 200
#foreach(sim=1:num_sims, .combine='cbind') %do% {
#    temp = df_output %>% arrange(desc(loglik))
#    temp[1, 'loglik'] = rnorm(1, temp[1, 'loglik'], temp[1, 'loglik_se'])
#    temp$lik = exp(temp$loglik - max(temp$loglik))
#    tot_sims = rmultinom(1, 1, temp$lik)
#    tot_sims
#} -> simulations_to_do
#df_output$num_sims = rowSums(as.data.frame(simulations_to_do))
#liksample <- df_output %>% filter(num_sims > 0) %>% 
#  select(beta1, beta2_1, beta2_2, beta2_3, num_init_1, num_init_2, num_init_3, loglik, loglik_se, lik, num_sims, nparticles, n_reps) 
#write.csv(liksample, file=paste0('../', liksample_filename))
#setwd('../')

design = read.csv(output_filename)

#### PREP FOR SIMULATION

simstart = convert_date_to_time(simstart, t_ref)
simend = convert_date_to_time(simend, t_ref)
min_data_time = convert_date_to_time(min_data_time, t_ref)
intervention_start = convert_date_to_time(intervention_start, t_ref)
min_data_time_ICU = convert_date_to_time(min_data_time_ICU, t_ref)

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
data_1 = get_civis_data(data_filename_1)  %>% mutate(ObsDeaths=ObsDeaths, ObsICU=ObsICU, Region=1) %>% select(time, ObsDeaths, ObsICU, Region)
data_2 = get_civis_data(data_filename_2)  %>% mutate(ObsDeaths=ObsDeaths, ObsICU=ObsICU, Region=2) %>% select(time, ObsDeaths, ObsICU, Region)
data_3 = get_civis_data(data_filename_3)  %>% mutate(ObsDeaths=ObsDeaths, ObsICU=ObsICU, Region=3) %>% select(time, ObsDeaths, ObsICU, Region)
data = as.data.frame(rbind(data_1, data_2, data_3)) %>%
  gather(Compartment, Cases, c(ObsDeaths, ObsICU)) %>%
  mutate(Time=time + as.Date('2020-01-14'))

print('Initializing parameters')
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)

pars$beta1 = design[1, 'beta1']
pars$beta2_1 = design[1, 'beta2_1']
pars$beta2_2 = design[1, 'beta2_2']
pars$beta2_3 = design[1, 'beta2_3']
pars$num_init_1 = round(design[1, 'num_init_1'])
pars$num_init_2 = round(design[1, 'num_init_2'])
pars$num_init_3 = round(design[1, 'num_init_3'])
pars$tmin=simstart
pars$tmax=simend

# add in invervention scaling
pars = add_interventions(covid_get_path(intervention_file), pars)
print(pars)

## Make a pomp object for inference
print('Simulating')

simout <- simulate_pomp_covid(
    n_regions = pars$n_regions,
    n_age_groups = n_age_groups,
    nsim=100,
    input_params = pars,
    delta_t = deltaT,
    contacts=pomp_contacts,
    population=population_list,
    nu_scales=nu_scales,
    beta_scales=beta_scales,
    frac_underreported=fraction_underreported,
    rprocess_Csnippet = rprocess_snippet,
    rinit_Csnippet = rinit_snippet
    ) %>% process_pomp_covid_output(agg_regions=F)

sims <- simout$plotting_output %>% ungroup() %>% filter(Time >= min_data_time)

get_obsprob = function(Time){
    pars$nu_3 * pars$theta_test * (1 - fraction_underreported[as.character(Time), 'frac_underreported'] * runif(length(Time), 0.8, 1))
}


# Impose observation model on hospitalized deaths
hd = sims %>% filter(Compartment == 'nHD') %>%
    mutate(Cases = rbetabinom(length(Time), Cases, get_obsprob(Time), pars$dispersion),
           Time = as.Date('2020-01-14') + Time)

ggplot() + 
  geom_line(hd, mapping=aes(x=Time, y=Cases, group=SimID), alpha=0.1) + 
  geom_point(data %>% filter(Compartment=='ObsDeaths'), mapping=aes(x=Time, y=Cases), color='red') +
  ylab('Observed hospitalized deaths') +
  facet_wrap(~Region, scales='free_y')

ggsave('fit_to_deaths.png', width=10, height=5)

icu = sims %>% filter(Compartment == 'IC') %>%
    mutate(Cases = rbetabinom(length(Time), Cases, get_obsprob(Time), pars$dispersion),
      Time = as.Date('2020-01-14') + Time)

ggplot() + 
  geom_line(icu, mapping=aes(x=Time, y=Cases, group=SimID), alpha=0.1) + 
  geom_point(data %>% filter(Compartment=='ObsICU'), mapping=aes(x=Time, y=Cases), color='red') +
  geom_vline(xintercept=as.Date('2020-04-07')) +
  ylab('Observed ICU cases') +
  facet_wrap(~Region, scales='free_y')
ggsave('fit_to_ICU.png', width=10, height=5)