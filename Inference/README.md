# Inference for Illinois SARS-CoV-2 model

We fitted the number of people initially infected in each region of Illinois (see [Data](../Data)) and the pre- and post-intervention transmission rates to the number of observed hospitalized deaths per day in each region.
Other parameters were fixed ([Parameters](.../Parameters)).

From a set of parameters, we can simulate both the hidden (latent) states (e.g., underlying fraction infected and immune) and observations (e.g., recorded COVID-19 deaths) of the model.
The likelihood of a parameter set measures how well it can reproduce the observed data.
To calculate the likelihood, we use a [particle filter](https://kingaa.github.io/sbied/pfilter/pfilter.html) implemented in the [`pomp` R package](https://kingaa.github.io/pomp/).
To find the best-fitting set of parameters, we use the [iterated filtering method](https://kingaa.github.io/sbied/mif/mif.html) also implemented in `pomp`.
Briefly, through repeated rounds of simulating the dynamics and then perturbing parameter values, this approach can find the parameter values that can best reproduce the data.

## Model file descriptions

* `rprocess_noise_non_hosp_deaths_region_contacts.c`: Process model that describes how people move through compartments at each timestep.
* `initializer_non_hosp_deaths.c`: Initializer that places people into infectious classes at the beginning of the simulation.
* `dmeasure_deaths_aggregate.c`: Calculation of model likelihood based on observed hospitalized deaths. Observed deaths at each timestep are a sample of those who have died in the hospital. We assume that we do not observe all COVID-19 deaths because testing does not detect all infections.
* `dmeasure_deaths_ICU_aggregate.c`: Calculation of model likelihood based on observed hospitalized deaths and confirmed ICU cases. Observed deaths at each timestep are a sample of those who have died in the hospital. Observed ICU cases are a sample of all people in ICU model compartments. We assume that we do not observe all COVID-19 deaths and ICU cases because testing does not detect all infections.

## Run script descriptions
* `inference_functions.R`: Essential functions for creating and using `pomp` objects. 
* `1_mif_single.R`: Runs `mif` search from Latin hypercube sample of points.
* `2_get_end_points_from_mif.R`: Aggregate endpoints of mif chains.
* `3_pfilter_end_points.R`: Evaluates likelihood on end `mif` points with `pfilter`.
* `4_pfilter_grid_search.R`: Performs a pfilter search on a Latin hypercube sample of points that surround the MLE.
* `5_simulate_and_plot.R`: Aggregates all points and simulates from MLE.
