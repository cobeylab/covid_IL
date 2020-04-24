# Inference for Illinois SARS-CoV-2 model

We infer key transmission parameters by fitting the model to data, with other model parameters fixed at values informed by existing literature (see [Parameters](../Parameters)). Specifically, our goal is to infer the number of people initially infected in each region of Illinois (see [Data](../Data)) along with the pre- and post-intervention transmission rates. We do this by finding the set of parameters that maximizes the likelihood of our data, the number of observed hospitalized deaths per day in each region.

Given a set of parameters, we can simulate both the hidden and observed states of the model. The likelihood of a given parameter set is then a measure of how well a given set of parameters can reproduce the observed data. To do this, we use a [particle filtering approach](https://kingaa.github.io/sbied/pfilter/pfilter.html) implemented in the [`pomp` R package](https://kingaa.github.io/pomp/). In order to find the best-fitting set of parameters, we use the [iterated filtering method](https://kingaa.github.io/sbied/mif/mif.html) also implemented in `pomp`. Briefly, through an iterated process of simulating the data and then perturbing parameter values, the best set of parameter values can be found.

## Model file descriptions

* `rprocess_interventionbeta_IH4.c`: Process model that describes how people move through compartments at each timestep.
* `initializer_compartment_distribute_IH4.c`: Initializer that places people into infectious classes at the beginning of the simulation.
* `dmeasure_deaths_aggregate.c`: Calculation of model likelihood based on observed hospitalized deaths. Observed deaths at each timestep are a sample of those who have died in the hospital. We assume that we do not observe all COVID-19 deaths because testing does not detect all infections.

## Run script descriptions
* `inference_functions.R`: Essential functions for creating and using `pomp` objects. 
* `1_pfilter_search.R`: Evaluates likelihood on a sample of points with pfilter.
* `2_mif_single.R`: Runs `mif` searches starting at the 200 best points from the initial `pfilter` search.
* `3_pfilter_mif_results.R`: Evaluates likelihood on output of `mif` searches.
