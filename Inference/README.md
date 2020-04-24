# Inference approach

While we use existing literature to inform the values of most of our model's parameters, some key parameters must be inferred by fitting our model to data. Specifically, our goal is to infer the number of people initially infected in each region along with the pre- and post-intervention transmission rates. We do this by finding the set of parameters that maximizes the likelihood of our data, the number of observed hospitalized deaths per day in each region.

Given a set of parameters, we can simulate both the hidden and observed states of the model. The likelihood of a given parameter set is then a measure of how well a given set of parameters can reproduce the observed data. To do this, we use a [particle filtering approach](https://kingaa.github.io/sbied/pfilter/pfilter.html) implemented in `pomp`. In order to find the best-fitting set of parameters, we use the [iterated filtering method](https://kingaa.github.io/sbied/mif/mif.html) also implemented in `pomp`. Briefly, through an iterated process of simulating the data and then perturbing parameter values, the best set of parameter values can be found.

# File descriptions

* `rprocess_interventionbeta_IH4.c`: Process model that describes how people move through compartments at each timestep.
* `initializer_compartment_distribute_IH4.c`: Initializer that places people into infectious classes at the beginning of the simulation.
* `dmeasure_deaths_aggregate.c`: Calculation of model likelihood based on observed hospitalized deaths. Observed deaths at each timestep are a sample of those who have died in the hospital. We assume that we do not observe all covid deaths because testing does not detect all infections.
