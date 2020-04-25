# Forecasting SARS-CoV-2 dynamics for the state of Illinois

Contributors (alphabetically listed): Phil Arevalo, Ed Baskerville, Spencer Carran, Sarah Cobey, Katelyn Gostic, Lauren McGough, Sylvia Ranjeva, & Frank Wen

Correspondence to: cobey@uchicago.edu, [Cobey lab](https://cobeylab.uchicago.edu)


### Model overview

This mathematical model infers past SARS-CoV-2 transmission rates in Illinois and can be used to forecast community spread, hospital and ICU burden, and mortality under current and hypothetical public health interventions.

This is a compartmental SEIR model.
Compartments consist of individuals who are susceptible (S); exposed and infected (but not yet infectious) (E); infectious, i.e., able to infect others (I); and recovered and immune (R).
We subdivide these compartments to track asymptomatic and symptomatic infections, multiple stages of hospitalization, and fatalities.
The model is age-structured, in that compartments are further subdivided into age groups that differ in, e.g., their probability of being an asymptomatic case and contacting other age groups.

When individuals are infected, they enter the latent or exposed class (E), where they cannot transmit infections to others (Figure 1).
Individuals that enter the asymptomatic class (A) can infect others in the community, but they eventually recover without showing symptoms.
People who enter the presymptomatic infectious class (P) eventually develop symptoms.
Symptomatic infections are divided into mild cases, I<sub>M</sub>, which resolve without hospital attention, and severe cases, I<sub>S</sub>, which require hospitalization.
Severe cases can be in the ICU (I<sub>C</sub) or not (class I<sub>H</sub>), and can die.
Deaths that occur outside of the hospital are estimated using a scaling factor informed by the model parameters and other information on the frequency of COVID-19 deaths outside hospitals (see [Parameters](./Parameters)).

![Figure 1](model_diagram.png)
Rate parameters determine the average amount of time in each compartment.
Probabilities determine the fraction of people following specific paths between compartments.
The model is stochastic, in that in each time step, the number of individuals transitioning between compartments is drawn randomly based on these rate and probability parameters.

### Data
The model is fitted to deaths reported by the New York Times from March 15 to March 24 and to in-hospital deaths reported by the Illinois Department of Public Health (IDPH) after March 24.
Deaths before March 15 are excluded due to concerns about excessive underreporting.
The New York Times data draw from publicly available data shared by IDPH.
The data we received from IDPH track in-hospital deaths precisely from March 24 onward.
To better approximate dynamics for the entire state, epidemic dynamics are estimated separately for three geographic subregions distinguished by the similarity of their epidemic activity to date ([Data](./Data)).

### Observation model
Limited testing capacity and false negatives mean that not all deaths from COVID-19 infection are observed.
Although the model tracks all underlying infections and deaths, it assumes only a fraction will be confirmed and counted.

### Inference
For each region in Illinois, we infer the transmission rate of SARS-CoV-2 before and during shelter-in-place.
The model also estimates the number of individuals infected on March 1, the start of the simulation.
Other parameters are fixed based on values from the literature ([Parameters](./Parameters)).
The model is fitted to the data using sequential Monte Carlo, a particle filter ([Inference](./Inference)).

### Model outputs
We simulate the dynamics of SARS-CoV-2 in Illinois using the best-fit parameters from the model inference.
The model can forecast different public health interventions, which will be later uploaded to the [Forecasting](./Forecasting) directory.
The forecasts incorporate several types of uncertainty that contribute to variation:
* Uncertainty in precisely how many individuals will be infected, recover, or die each day (demographic stochasticity)
* Uncertainty in the fraction of COVID-19 deaths in hospitals that were missed in March due to inadequate testing
* Uncertainty in the fraction of the population infected on March 1
* Uncertainty in the inferred transmission rates
Incorporating all of these types of uncertainty results in a wide range of potential epidemic trajectories.

### Public health interventions 
The model can be adapted to investigate different types of interventions.
We currently model shelter-in-place as a reduction in the overall transmission rate, which is shared by all infected individuals, and reductions in age- and location-specific contact rates.
The three scenarios we explore are:

1. Shelter in place extended indefinitely ("indefinite"). The reduction in the overall transmission rate is inferred (estimated). We assume a 100% reduction in at-school contacts, a 40% reduction in at-work contacts, and a 50% reduction in contacts outside of the home.
2. Shelter in place is lifted on May 1, 2020 ("lifted")
3. Shelter in place is never adopted, and individuals do not change their behavior from early March ("never")

### Caveats

The model's predictions will shift as new data on COVID-19 emerges from Illinois and around the world.
We will update the model to incorporate better assumptions about the underlying biology and epidemiology.
We will also continue to try different modeling approaches with different assumptions to see how they affect conclusions.
And we will continue to engage in a dialogue with other scientists and modelers, and learn from what they are doing (e.g., via the [MIDAS Network](https://midasnetwork.us/) and the [Mobility Data Network](https://www.covid19mobility.org/)).

Some important assumptions of our model include:
* Individuals are effectively no longer infectious once they enter hospitals. Obviously, if insufficient PPE is available or a case is not diagnosed in time, this could be a bad assumption.
* The geographic regions of Illinois are independent (i.e., meaningful transmission does not occur between regions).
* Interventions occur immediately at a fixed intensity. This assumption is currently supported by mobility data we have analyzed from Facebook.
* Many parameters in our model are fixed based on existing literature (see [Parameters](./Parameters)).

### Reference
1. King AA, Nguyen D and Ionides EL (2015) Statistical inference for partially observed Markov processes via the R package pomp. arXiv preprint arXiv:1509.00503.

### License for text and figures

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />All text and figures in this repository are licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.

### License for code

All code in this repository is Copyright © 2020.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
