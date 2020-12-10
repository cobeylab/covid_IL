simend = '2020-12-07' # Set to last date in the EMR data
project_end = '2021-04-01'
project_zoom = as.Date('2020-03-01')
model_name = 'beta_hfr_mu_gamma_constphi_cdc'

## Data to fit
ll_file = 'Inference/projections_20201209/data/ll_deaths_20201208.csv'
emr_data_file = 'Inference/projections_20201209/data/emr_fit_long_20201208.csv'
emr_report = 'Inference/projections_20201209/data/emr_report_long_20201208.csv'
hosp_capacity_file = 'Inference/projections_20201209/data/capacity_weekday_average_20201208.csv'
debug.mif = F
start_fit_at_mle = T
mle_file =sprintf('old_mle_agg_mif_%s.csv', model_name)

## Parameters for global search and likelihood calculation, 
## set all of these to 2 if you want to just test that the pipeline works
## Otherwise, n_mif = 75, n_particles_mif = 3000, and n_particles_pfilter = 5000
n_mif = 25
n_particles_mif = 8000
n_particles_pfilter = 8000
n_reps_pfilter=5
cooling_rate = 0.99
sd_val = 0.01
n_sim = 1000
if (debug.mif == T){
    n_mif = 2
    n_particles_mif = 2
    n_particles_pfilter = 2
    n_reps_pfilter=2
    n_sim = 2
}

use_mle = T
output_dir = sprintf('./inference_output_%s/',model_name)

# CSV that contains columns for the time and each of the covariates that you want to use
# Note that the filepath should be specified RELATIVE to the covid_IL directory.
beta_covariate = 'Inference/projections_20201209/data/beta_values_no_age.covar.cli_tpr.csv'
beta_covariate_column = 'None' # Name of the transmission covariate you want to use
beta_scale = 'Data/beta_scaling_mobility.csv'
default_par_file = './default_parameter_values.csv'

### Shouldn't need to change anything below this point ###

t_ref='2020-01-14'
simstart = '2020-03-01'
min_data_time = '2020-03-16'
intervention_start='2020-03-16'
min_data_time_ICU = '2020-04-03'

idph_filename = './idph_public_covid_region.csv' # Location of IDPH public data broken down by region
dmeasFile = 'Inference/projections_20201209/dmeasure_no_age.c'
init_file='Inference/projections_20201209/initializer_no_age.c'
rmeasFile = 'Inference/projections_20201209/rmeasure_no_age.c'
rprocFile = 'Inference/projections_20201209/rprocess_no_age_changepoint.c'
function_file='Inference/projections_20201209/simulation_functions_no_age.R'

# all covariates
fraction_underreported_file='internal_data/2020-11-23_IL_underrerporting.csv'
nonhosp_deaths = 'internal_data/2020-10-24_IL_frac_deaths_nonhosp.csv'
population_file = 'Data/covid_region_populations.csv'
cli_file = 'internal_data/cli_admissions_2020-11-24.csv'

print('finished reading input files')