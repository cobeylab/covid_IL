## Specify model input files
initFile = '../Inference/initializer_non_hosp_deaths.c'
rprocFile = '../Inference/rprocess_calculate_phis.c'
rmeasFile = '../Inference/rmeasure_hosp_non_hosp_ICU.c'
nu_scales_file = '../Data/nu_scaling.csv'
fraction_underreported_file = "../Data/frac_underreported.csv"  
default_par_file = "./evaluate_model_changes/default_parameter_values_scenario6.csv"
contact_filename = '../Data/IL_4metaregions_symmetric.RData'
intervention_table_filename = "./intervention_table.csv"

## Specify population filenames
population_filename_1= "../Data/IL_population_north-central.csv"
population_filename_2= "../Data/IL_population_central.csv"
population_filename_3= "../Data/IL_population_northeast.csv"
population_filename_4= "../Data/IL_population_southern.csv"

IDPH_death_data_filename = "../Data/covidtracking_IDPH.csv" # Must include file with IDPH reported new deaths per day
regional_idph_file = "../Data/incident_idph_regions.csv"