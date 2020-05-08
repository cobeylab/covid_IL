## Specify model input files
initFile = '../Inference/initializer_non_hosp_deaths.c'
rprocFile = '../Inference/rprocess_noise_non_hosp_deaths_region_contacts.c'
nu_scales_file = '../Data/nu_scaling.csv'
fraction_underreported_file = "../Data/frac_underreported.csv"  
default_par_file = "../Parameters/parameter_values.csv"
contact_filename = '../Data/IL_POMP_contacts_expanded.RData'
intervention_table_filename = "intervention_table.csv"

## Specify population filenames
population_filename_1= "../Data/IL_population_EMS_regions_1-2.csv"
population_filename_2= "../Data/IL_population_EMS_regions_3-6.csv"
population_filename_3= "../Data/IL_population_EMS_regions_7-11.csv"
