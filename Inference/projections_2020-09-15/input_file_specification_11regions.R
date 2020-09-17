simend = '2020-09-14' # Set to last date in the EMR data
use_changepoint  = T # Model changes in transmission as changepoints or with covariates?
beta_covariate_column = 'None' # Name of the transmission covariate you want to use

## Parameters for global search and likelihood calculation, 
## set all of these to 2 if you want to just test that the pipeline works
## Otherwise, n_mif = 75, n_particles_mif = 3000, and n_particles_pfilter = 5000
n_mif = 50
n_particles_mif = 3000
n_particles_pfilter = 5000

# CSV that contains columns for the time and each of the covariates that you want to use
# Note that the filepath should be specified RELATIVE to the covid_IL directory.
beta_covariate_file = 'Data/mobility_covar_table.csv' 

### Shouldn't need to change anything below this point ###

t_ref='2020-01-14'
simstart = '2020-03-01'
min_data_time = '2020-03-16'
intervention_start='2020-03-16'
min_data_time_ICU = '2020-04-07'

covid_region_data_filename='/project2/cobey/covid-civis/processed_data/emr_linelist_fitting_data_covid.csv'
idph_filename = '/project2/cobey/covid-civis/processed_data/idph_public/idph_public_covid_region.csv' # Location of IDPH public data broken down by region
dmeasFile = 'Inference/projections_2020-09-15/dmeasure_project11.c'
init_file='Inference/projections_2020-09-15/initializer_project11.c'
if (use_changepoint){
    rprocFile='Inference/projections_2020-09-15/rprocess_project11.c' 
}
rmeasFile = 'Inference//projections_2020-09-15/rmeasure_project11.c'
inference_file='Inference/projections_2020-09-15/inference_functions.R'

fraction_underreported_file='Data/2020-08-11_IL_underrerporting.csv'
contact_filename='Data/contacts_covar_table.csv'
intervention_file='Forecasting/intervention_table.csv'
age_dist_file = 'Parameters/consolidated_age_distributions.csv'

## Projections
projection_dir = paste0(Sys.Date(), "/")
final_projection_date = "2021-03-31"
projection_functions = "Forecasting/projection_functions.R"
model_name="baseline"
population11_file = 'Data/covid_region_populations.csv'
