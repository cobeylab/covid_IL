simend = '2020-07-22' # Set to last date in the EMR data
use_changepoint  = F # Model changes in transmission as changepoints or with covariates?
beta_covariate_column = 'crowdedness' # Name of the transmission covariate you want to use

## Parameters for global search and likelihood calculation, 
## set all of these to 2 if you want to just test that the pipeline works
## Otherwise, n_mif = 75, n_particles_mif = 3000, and n_particles_pfilter = 5000
n_mif = 2
n_particles_mif = 2
n_particles_pfilter = 2

# CSV that contains columns for the time and each of the covariates that you want to use
# Note that the filepath should be specified RELATIVE to the covid_IL directory.
beta_covariate_file = 'Data/mobility_covar_table.csv' 

### Shouldn't need to change anything below this point ###

t_ref='2020-01-14'
simstart = '2020-03-01'
min_data_time = '2020-03-16'
intervention_start='2020-03-16'
min_data_time_ICU = '2020-04-07'

data_filename = './emr_linelist_fitting_data.csv' # Location of linelist data
idph_filename = './incident_idph_regions.csv' # Location of IDPH public data broken down by region
dmeasFile = 'Inference/dmeasure.c'
init_file='Inference/initializer.c'
if (use_changepoint){
    rprocFile='Inference/rprocess_changepoint.c' 
} else{
    rprocFile='Inference/rprocess_covariate_beta.c'
}

rmeasFile = 'Inference/rmeasure.c'

fraction_underreported_file='Data/frac_underreported.csv'
contact_filename='Data/contacts_covar_table.csv'
inference_file='Inference/inference_functions.R'
intervention_file='Forecasting/intervention_table.csv'
age_dist_file = 'Parameters/consolidated_age_distributions.csv'

population_filename_1 = 'Data/IL_population_north-central.csv'
population_filename_2 = 'Data/IL_population_central.csv'
population_filename_3 = 'Data/IL_population_northeast.csv'
population_filename_4 = 'Data/IL_population_southern.csv'
population_filename_5 = 'Data/IL_population_chicago.csv'

## Projections
projection_dir = paste0('./projections_', Sys.Date(), "/")
final_projection_date = "2021-03-13"
projection_functions = "Forecasting/projection_functions.R"
model_name="baseline"