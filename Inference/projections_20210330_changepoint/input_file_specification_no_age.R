
simend = '2021-03-28' # Set to last date in the EMR data
project_end = '2021-09-01' # Last date for projections
min_data_time_ICU = '2020-04-08' # Time when ICU data starts
t_ref='2020-01-14' # Defined as t = 1
simstart = '2020-03-01' # When to start pomp
min_data_time = '2020-03-16' # When do the death data start
intervention_start='2020-03-16' # Deprecated
project_zoom = as.Date('2021-01-01') # minimum time for zoomed in view

model_name = 'changepoints'
previous_inference = '../projections_20210316_changepoint/inference_output_changepoints/'
homedir = 'Inference/projections_20210330_changepoint/'

debug.mif = F
start_fit_at_mle = T
simulate_new = T
use_mle = T

## Data to fit
ll_file = sprintf('internal_data/ll_deaths_2021-03-30.csv')
emr_data_file = sprintf('internal_data/emr_fit_long_2021-03-30.csv')
emr_report = sprintf('internal_data/emr_report_long_2021-03-30.csv')
hosp_capacity_file = sprintf('internal_data/capacity_weekday_average_20210330.csv')
hospital_covar_file = sprintf('internal_data/hospital_covariates_2021-03-02.csv')
frac_averted_filename = 'internal_data/frac_averted_2021-03-30.csv'
cli_file = 'internal_data/cli_admissions_2021-03-25.csv'
fraction_underreported_file='internal_data/2020-01-17_IL_underrerporting.csv'

population_file = 'Data/covid_region_populations.csv'

mle_file =sprintf('old_mle_agg_mif_%s.csv', model_name)
new_mle_file=sprintf('agg_mif_%s.csv', model_name)

n_mif = 30
n_particles_mif = 8000
n_particles_pfilter = 8000
n_reps_pfilter=5
cooling_rate = 0.99
sd_val = 0.02
sd_val_init = 0.02
n_sim = 1000
if (debug.mif == T){
    n_mif = 2
    n_particles_mif = 2
    n_particles_pfilter = 2
    n_reps_pfilter=5
    n_sim = 2
}

output_dir = sprintf('./inference_output_%s/', model_name)
modeldir = sprintf('Inference/model_files/%s/', model_name)

beta_covariate = ''
beta_covariate_column = 'None' # Name of the transmission covariate you want to use
default_par_file = './default_parameter_values.csv'

dmeasFile = sprintf('%s/dmeasure.c', modeldir)
init_file= sprintf('%s/initializer.c', modeldir)
rmeasFile = sprintf('%s/rmeasure.c', modeldir)
rprocFile = sprintf('%s/rprocess.c', modeldir)
rprocScenarioFile = sprintf('%s/rprocess_scenario.c', modeldir)
globalFile = sprintf('%s/globals.c', modeldir)
function_file=sprintf('%s/functions.R', modeldir)

print('finished reading input files')

library(zoo)