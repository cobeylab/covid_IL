# Inference for Illinois SARS-CoV-2 model

## Summary

We calibrated our model to infer the number of people initially infected in each region of Illinois (see [Forecasts](../Forecasts)); region-specific transmission rates at 13 changepoints; the region-specific, time-varying hospitalization fatality rate (HFR); the region-specific, time-varying duration of the hospital stay; and the region-specific fraction of non-hospitalized infections that result in death.
Other parameters were fixed ([Parameters](.../Parameters)).

From a set of parameters, we can simulate both the hidden (latent) states (e.g., underlying fraction infected and immune) and observations (e.g., recorded COVID-19 deaths) of the model.
The likelihood of a parameter set measures how well it can reproduce the observed data.
To calculate the likelihood, we use a [particle filter](https://kingaa.github.io/sbied/pfilter/pfilter.html) implemented in the [`pomp` R package](https://kingaa.github.io/pomp/).
To find the best-fitting set of parameters, we use the [iterated filtering method](https://kingaa.github.io/sbied/mif/mif.html) also implemented in `pomp`.
Briefly, through repeated rounds of simulating the dynamics and then perturbing parameter values, this approach can find the parameter values that can best reproduce the data.

Code used to fit the model to new data are stored in dated directories. The most recent fits were produced with the code in [projections_20201124](./projections_20201124).

## Model file descriptions

* `rprocess_no_age_changepoint.c`: Process model that describes how people move through compartments at each timestep.
* `initializer_no_age.c`: Initializer that places people into infectious classes at the beginning of the simulation.
* `dmeasure_no_age.c`: Calculation of model likelihood based on observed hospitalized deaths, total observed deaths, and census hospitalization.
* `rmeasure_no_age.c`: Code to simulate observed states from latent states.

## Run script descriptions
* `simulation_functions_no_age.R`: Essential functions for creating and using `pomp` objects. 
* `fit_regions.R`: Runs `mif` search starting from previous MLE.
* `aggregate_results.R`: Aggregate endpoints of mif chains.
* `simulate_regions.R`: Simulate results from MLE of each region.
* `slice_regions.R`: Take a likelihood slice over the final transmission rate.
* `project_regions.R`: Make final projections.
* `project_statewide.R`: Aggregate results from individual regions to generate Illinois projections.

