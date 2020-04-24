# Forecasting for Illinois COVID-19 model
## File descriptions

* `simulate_projections.Rmd`: R-markdown file to simulate and plot forecasts for the state of IL based on the inference results.
* `input_file_specification.R`: R script to read in files for simulation. Called within `simulate_projections.Rmd`
* `inference_to_simulation.R`: R script to run and store simulation output. Called within `simulate_projections.Rmd`
* `simulation_statewide.R`: R script with functions for simulation. Called within `inference_to_simulation.R`
* `plot_intervention_comparisons.R`: R script to plot simulation output. Called within `simulate_projections.Rmd`
