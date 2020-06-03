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

root <- '../../'
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
data_filename = args[3]
population_filename_1 = args[4]
population_filename_2 = args[5]
population_filename_3 = args[6]
simstart = args[7]
simend = args[8]
min_data_time = args[9]
intervention_start= args[10]
num_cores=as.numeric(args[11])
maxjobs=as.numeric(args[12])
min_data_time_ICU = args[13]
init_file=args[14]
rprocFile=args[15]
nu_scales_file=args[16]
fraction_underreported_file=args[17]
contact_filename=args[18]
inference_file=args[19]
t_ref=args[20]
intervention_file=args[21]
population_filename_4 = args[22]
region_to_test=as.numeric(args[23])

source(covid_get_path(inference_file))

# Aggregate points
print('Aggregating points')
output_filename <- paste0(output_dir, "consolidated.all.pfilter.csv")
liksample_filename <- paste0(output_dir, "consolidated.all.sample_from_surface.csv")
file_path = output_dir #Replace this with a file path to the directory storing the csv files 
setwd(file_path)
file_list <- list.files(pattern = "output_pfilter_.*.csv$")
data <- rbindlist(lapply(file_list, fread))
df_output <- as.data.frame(data) %>% arrange(desc(loglik)) %>%
    select(beta1_1,beta1_2,beta1_3,beta1_4, beta2_1, beta2_2, beta2_3, beta2_4, region_non_hosp_1, region_non_hosp_2, region_non_hosp_3, region_non_hosp_4, num_init_1, num_init_2, num_init_3, num_init_4, b_elderly, loglik, loglik_se, nparticles, n_reps) 
write.csv(df_output, file = paste0('../', output_filename))

df_output$lik = exp(df_output$loglik - max(df_output$loglik))
num_sims = 200
foreach(sim=1:num_sims, .combine='cbind') %do% {
    temp = df_output %>% arrange(desc(loglik))
    temp[1, 'loglik'] = rnorm(1, temp[1, 'loglik'], temp[1, 'loglik_se'])
    temp$lik = exp(temp$loglik - max(temp$loglik))
    tot_sims = rmultinom(1, 1, temp$lik)
    tot_sims
} -> simulations_to_do
df_output$num_sims = rowSums(as.data.frame(simulations_to_do))
liksample <- df_output %>% filter(num_sims > 0) %>% 
  select(beta1_1,beta1_2,beta1_3,beta1_4, beta2_1, beta2_2, beta2_3, beta2_4, region_non_hosp_1, region_non_hosp_2, region_non_hosp_3, region_non_hosp_4, num_init_1, num_init_2, num_init_3, num_init_4, b_elderly, loglik, loglik_se, lik, num_sims, nparticles, n_reps) 
write.csv(liksample, file=paste0('../', liksample_filename))
setwd('../')

design = read.csv(output_filename)
print(design)
source('set_up_covariates_and_data.R')

pars$beta1_1 = design[1, 'beta1_1']
pars$beta1_2 = design[1, 'beta1_2']
pars$beta1_3 = design[1, 'beta1_3']
pars$beta1_4 = design[1, 'beta1_4']

pars$beta2_1 = design[1, 'beta2_1']
pars$beta2_2 = design[1, 'beta2_2']
pars$beta2_3 = design[1, 'beta2_3']
pars$beta2_4 = design[1, 'beta2_4']

pars$region_non_hosp_1 = design[1, 'region_non_hosp_1']
pars$region_non_hosp_2 = design[1, 'region_non_hosp_2']
pars$region_non_hosp_3 = design[1, 'region_non_hosp_3']
pars$region_non_hosp_4 = design[1, 'region_non_hosp_4']

pars$num_init_1 = round(design[1, 'num_init_1'])
pars$num_init_2 = round(design[1, 'num_init_2'])
pars$num_init_3 = round(design[1, 'num_init_3'])
pars$num_init_4 = round(design[1, 'num_init_4'])

pars$b_elderly = design[1, 'b_elderly']

pars$tmin=simstart
pars$tmax=simend
print(pars)

## Make a pomp object for inference
print('Simulating')

simout <- simulate_pomp_covid(
    n_regions = pars$n_regions,
    n_age_groups = n_age_groups,
    nsim=50,
    input_params = pars,
    delta_t = deltaT,
    contacts=pomp_contacts,
    population=population_list,
    nu_scales=nu_scales,
    beta_scales=beta_scales,
    frac_underreported=fraction_underreported,
    rprocess_Csnippet = rprocess_snippet,
    rinit_Csnippet = rinit_snippet,
    rmeasure_Csnippet=rmeasure_snippet,
    obsnames=c(paste0('ObsHospDeaths_', 1:4), paste0('ObsNonHospDeaths_', 1:4), paste0('ObsICU_', 1:4)) 
    ) %>% process_pomp_covid_output(agg_regions=F)

regions = c('northcentral', 'central', 'northeast', 'southern')
sims <- simout$plotting_output %>% ungroup()
#print(head(sims))
sims$Region = sapply(sims$Region, FUN=function(x){regions[as.numeric(x)]})
#print(head(sims))

data = read.csv(data_filename)
names(data) = c('date', 'Region', 'ObsHospDeaths', 'ObsNonHospDeaths', 'ObsICU', 'Source')
data = data %>% gather(Compartment, Cases, c(ObsHospDeaths, ObsNonHospDeaths, ObsICU)) %>%
  mutate(Time=as.Date(date)) %>% filter(Region %in% regions)

# Impose observation model on hospitalized deaths
hd = sims %>% filter(Compartment == 'Reported hd') %>%
    mutate(Time = as.Date('2020-01-14') + Time)

ggplot() + 
  geom_line(hd, mapping=aes(x=Time, y=Cases, group=SimID), alpha=0.1) + 
  geom_point(data %>% filter(Compartment=='ObsHospDeaths'), mapping=aes(x=Time, y=Cases), color='red') +
  ylab('Observed hospitalized deaths') +
  facet_wrap(~Region, scales='free_y')

ggsave(paste0(output_dir, 'fit_to_hosp_deaths.png'), width=10, height=5)


nhd = sims %>% filter(Compartment == 'Reported nhd') %>%
    mutate(Time = as.Date('2020-01-14') + Time)

ggplot() + 
  geom_line(nhd, mapping=aes(x=Time, y=Cases, group=SimID), alpha=0.1) + 
  geom_point(data %>% filter(Compartment=='ObsNonHospDeaths'), mapping=aes(x=Time, y=Cases), color='red') +
  ylab('Observed non-hospitalized deaths') +
  facet_wrap(~Region, scales='free_y')

ggsave(paste0(output_dir, 'fit_to_nonhosp_deaths.png'), width=10, height=5)

icu = sims %>% filter(Compartment == 'ObsICU') %>%
    mutate(Time = as.Date('2020-01-14') + Time)

ggplot() + 
  geom_line(icu, mapping=aes(x=Time, y=Cases, group=SimID), alpha=0.1) + 
  geom_point(data %>% filter(Compartment=='ObsICU'), mapping=aes(x=Time, y=Cases), color='red') +
  geom_vline(xintercept=as.Date('2020-04-07')) +
  ylab('Observed ICU cases') +
  facet_wrap(~Region, scales='free_y')
ggsave(paste0(output_dir, 'fit_to_ICU.png'), width=10, height=5)