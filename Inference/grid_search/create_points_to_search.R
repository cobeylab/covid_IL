library(pomp)
library(tidyverse)

region = 3

# initialize design
num_points = 40000

# initialize design
lower_pars = c(beta1_1=0.01,
    beta1_2=0.01,
    beta1_3=0.02,
    beta1_4=0.01, 
  beta2_1=0.01, 
  beta2_2=0.01, 
  beta2_3=0.01,
  beta2_4=0.01,
  num_init_1=10,
  num_init_2=10,
  num_init_3=5000,
  num_init_4=10,
  region_non_hosp_1=0.2,
  region_non_hosp_2=0.2,
  region_non_hosp_3=0.2,
  region_non_hosp_4=0.2,
  b_elderly=0)

upper_pars = c(beta1_1=0.05,
    beta1_2=0.05,
    beta1_3=0.045,
    beta1_4=0.05, 
  beta2_1=0.03, 
  beta2_2=0.03, 
  beta2_3=0.025,
  beta2_4=0.03,
  num_init_1=800,
  num_init_2=800,
  num_init_3=20000,
  num_init_4=800,
  region_non_hosp_1=0.5,
  region_non_hosp_2=0.5,
  region_non_hosp_3=0.5,
  region_non_hosp_4=0.5,
  b_elderly=0.6)

if (region != -1){
  param_names = names(lower_pars)[grepl(paste0('_',region), names(lower_pars))]

  lower = lower_pars[param_names]
  upper = upper_pars[param_names]
  #lower_pars=lower_pars[c(param_names, 'b_elderly')]
  #upper_pars=upper_pars[c(param_names, 'b_elderly')]

  design = expand.grid(
    beta1=seq(lower[1], upper[1], length.out=10),
  beta2=seq(lower[2], upper[2], length.out=10),
  num_init=seq(lower[3], upper[3], length.out=40),
  region_non_hosp=seq(lower[4], upper[4], length.out=5),
  b_elderly=seq(0, 0.5, length.out=5))
  
  design = design %>% filter(beta1 > beta2) %>% mutate(region_to_test=region)
  names(design) = c(param_names, 'b_elderly')
}

#design = sobolDesign(lower=lower_pars, 
# upper=upper_pars, 
#  num_points)

write.csv(design, sprintf('grid_search_params.csv',region), row.names=F)