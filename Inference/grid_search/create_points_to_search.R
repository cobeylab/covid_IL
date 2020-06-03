library(pomp)

# initialize design
num_points = 40000

# initialize design
lower_pars = c(beta1_1=0.005,
    beta1_2=0.005,
    beta1_3=0.005,
    beta1_4=0.005, 
  beta2_1=0.005, 
  beta2_2=0.005, 
  beta2_3=0.005,
  beta2_4=0.005,
  num_init_1=10,
  num_init_2=10,
  num_init_3=100,
  num_init_4=10,
  region_non_hosp_1=0.2,
  region_non_hosp_2=0.2,
  region_non_hosp_3=0.2,
  region_non_hosp_4=0.2,
  b_elderly=0)

upper_pars = c(beta1_1=0.05,
    beta1_2=0.05,
    beta1_3=0.05,
    beta1_4=0.05, 
  beta2_1=0.05, 
  beta2_2=0.05, 
  beta2_3=0.05,
  beta2_4=0.05,
  num_init_1=1000,
  num_init_2=1000,
  num_init_3=15000,
  num_init_4=1000,
  region_non_hosp_1=0.7,
  region_non_hosp_2=0.7,
  region_non_hosp_3=0.7,
  region_non_hosp_4=0.7,
  b_elderly=1)

design = sobolDesign(lower=lower_pars, 
  upper=upper_pars, 
  num_points)

write.csv(design, 'points_to_search.csv', row.names=F)