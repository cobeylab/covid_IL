## Specify model input files
initFile = '../Inference/initializer.c'
rprocFile = '../Inference/rprocess.c'
rmeasFile = '../Inference/rmeasure.c'
nu_scales_file = '../Data/nu_scaling.csv'
fraction_underreported_file = "../Data/frac_underreported.csv"  
icu_reporting_file = "../Data/icu_reporting.csv"  
default_par_file = "./default_parameter_values.csv"
contact_filename = '../Data/IL_5metaregions_symmetric.RData'
intervention_table_filename = "./intervention_table.csv"

## Specify population filenames
population_filename_1= "../Data/IL_population_north-central.csv"
population_filename_2= "../Data/IL_population_central.csv"
population_filename_3= "../Data/IL_population_northeast.csv"
population_filename_4= "../Data/IL_population_southern.csv"
population_filename_5= "../Data/IL_population_chicago.csv"

IDPH_death_data_filename = "../Data/covidtracking_IDPH.csv" # Must include file with IDPH reported new deaths per day
regional_idph_file = "../Data/incident_idph_regions.csv"