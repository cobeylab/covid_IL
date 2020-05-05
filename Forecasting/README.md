# Forecasting for Illinois SARS-CoV-2 model

## Evaluating the effects of public health interventions
We implement public health interventions both by scaling contact rates and by scaling the transmission rate directly.
We assume that a stay-at-home order scales down non-home contacts (work contacts, school contacts, and other contacts). 
On March 16, we scaled down school contacts to 0, work contacts to 0.6 times their normal value, and other contacts to 0.5 times their normal value. 
We then inferred a separate post-intervention transmission rate.
For the baseline model scenario, we maintain the inferred value of the post-intervention transmission rate indefinitely. We assume that this is the maximum reduction in transmission. To model the relaxing of shelter-in-place interventions, we project three additional scenarios with a scaled increase to the post-intervention transmission rate by 20%, 40% or 60%.

## File descriptions

* `simulate_projections.Rmd`: R-markdown file to simulate and plot forecasts for the state of IL based on the inference results.
* `input_file_specification.R`: R script to read in files for simulation. Called within `simulate_projections.Rmd`
* `inference_to_simulation.R`: R script to run and store simulation output. Called within `simulate_projections.Rmd`
* `simulation_statewide.R`: R script with functions for simulation. Called within `inference_to_simulation.R`
* `plot_intervention_comparisons.R`: R script to plot simulation output. Called within `simulate_projections.Rmd`
