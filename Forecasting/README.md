# Forecasting for Illinois SARS-CoV-2 model

## Evaluating the effects of public health interventions
We implement public health interventions both by scaling contact rates and by scaling the transmission rate directly.
We assume that a stay-at-home order scales down non-home contacts (work contacts, school contacts, and other contacts). 
On March 16, we scaled down school contacts to 0, work contacts to 0.6 times their normal value, and other contacts to 0.5 times their normal value. 
We then inferred a separate post-intervention transmission rate.
For the baseline model scenario, we maintain the inferred value of the post-intervention transmission rate indefinitely. We assume that this is the maximum reduction in transmission. To model the relaxing of shelter-in-place interventions, we project three additional scenarios:

1. The post-intervention transmission rate increases by 20%. At our current parameter estimates, this corresopnds to a 72% reduction in pre-intervention transmission rate.
2. The post-intervention transmission rate increases by 40%. At our current parameter estimates, this corresopnds to a 67% reduction in pre-intervention transmission.
3. The post-intervention transmission rate increases by 60%. At our current estimates, this corresopnds to a 62% reduction in pre-intervention transmission.

Figure 1 gives the projected incident infections, population prevalence, non-ICU hospital occupancy, and ICU occupancy through June 13, 2020 under the four model scenarios. 

![Figure 1](./plots/summary_outputs.png)


Figure 2 gives the projected daily reported hosptialized and non-hospitalized deaths through June 13, 2020, and reported daily deaths from the Illinois Department of Public Health (IDPH) under the four model scenarios. 

![Figure 2](./plots/death_summary_outputs.png) 

## File descriptions

* `simulate_projections.Rmd`: R-markdown file to simulate and plot forecasts for the state of IL based on the inference results.
* `input_file_specification.R`: R script to read in files for simulation. Called within `simulate_projections.Rmd`
* `inference_to_simulation.R`: R script to run and store simulation output. Called within `simulate_projections.Rmd`
* `simulation_statewide.R`: R script with functions for simulation. Called within `inference_to_simulation.R`
* `plot_intervention_comparisons.R`: R script to plot simulation output. Called within `simulate_projections.Rmd`
