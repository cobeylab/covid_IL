data_filename = './emr_linelist_fitting_data.csv' # Location of linelist data
idph_filename = './incident_idph_regions.csv' # Location of IDPH public data broken down by region
simend = '2020-07-13' # Set to last date in data

### Shouldn't need to change anything below this poin ###
t_ref='2020-01-14'
simstart = '2020-03-01'
min_data_time = '2020-03-16'
intervention_start='2020-03-16'
min_data_time_ICU = '2020-04-07'

dmeasFile = 'Inference/dmeasure.c'
init_file='Inference/initializer.c'
rprocFile='Inference/rprocess.c'
rmeasFile = 'Inference/rmeasure.c'

fraction_underreported_file='Data/frac_underreported.csv'
contact_filename='Data/IL_5metaregions_symmetric.RData'
inference_file='Inference/inference_functions.R'
intervention_file='Forecasting/intervention_table.csv'
age_dist_file = 'Parameters/consolidated_age_distributions.csv'

population_filename_1 = 'Data/IL_population_north-central.csv'
population_filename_2 = 'Data/IL_population_central.csv'
population_filename_3 = 'Data/IL_population_northeast.csv'
population_filename_4 = 'Data/IL_population_southern.csv'
population_filename_5 = 'Data/IL_population_chicago.csv'