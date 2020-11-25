# Forecasting SARS-CoV-2 dynamics for the state of Illinois

Contributors from the [Cobey lab](https://cobeylab.uchicago.edu) (listed alphabetically): Phil Arevalo, Ed Baskerville, Spencer Carran, Sarah Cobey, Katelyn Gostic, Lauren McGough, Sylvia Ranjeva, & Frank Wen

### Model overview

This mathematical model infers past SARS-CoV-2 transmission rates in Illinois and can be used to forecast community spread, hospital burden, and mortality.

This is a compartmental SEIR model.
Compartments consist of individuals who are susceptible (S); exposed and infected (but not yet infectious) (E); infectious, i.e., able to infect others (I); and recovered and immune (R).
We subdivide these compartments to track infections, hospitalization, and fatalities.

When individuals are infected, they enter the latent or exposed class (E), where they cannot transmit infections to others (Figure).
People who enter the presymptomatic infectious class (P) either develop symptoms or remain asymptomatic.
We assume that the duration of infectiousness for asymptomatic infections is the same as the duration of infectiousness for mild symptomatic infections, and therefore do not explicitly model an asymptomatic compartment.
Presymptomatic infections are divided into non-hospitalized infections (I<sub>R</sub> and I<sub>D</sub>) and infections that require hospitalization (I<sub>H</sub>).
Hospitalized infections are either discharged (H<sub>R</sub>) or die in the hospital (H<sub>D</sub>).

![Figure 1](model_diagram.png)
Rate parameters (blue) determine the average amount of time in each compartment.
Probabilities (red) determine the fraction of people following specific paths between compartments.
The model is stochastic, in that in each time step, the number of individuals transitioning between compartments is drawn randomly based on these rate and probability parameters.

### Data
The model is fitted to in-hospital deaths reported by the Illinois Department of Public Health (IDPH) after March 15, all COVID deaths reported by IDPH after March 15, and the total number of hospital beds occupied by confirmed COVID cases reported by IDPH after May 6.
The data we received from IDPH track in-hospital deaths precisely from March 16 onward and confirmed cases in the hospital from May 7 onward.
To better approximate dynamics for the entire state, epidemic dynamics are estimated separately for 11 geographic subregions used by the [Restore Illinois plan](https://coronavirus.illinois.gov/s/restore-illinois-introduction) ([Data](./Data)).

### Observation model
Incomplete hospital reporting means that not all deaths and hospital admissions from COVID-19 infection are observed.
Although the model tracks all underlying infections and deaths, it assumes only a fraction will be confirmed and counted.

### Inference
For each region in Illinois, we infer the transmission rate at 13 timepoints and assume that the transmission rate changes linearly between each timepoint.
We also infer region-specific parameters governing the hospital course: the infection hospitalization ratio (IHR), the time-varying hospitalization fatality ratio (HFR) and the time-varying duration of the hospital stay.
To account for deaths that occur without hospitalization, we infer a region-specific probability that a non-hospitalized infection ends in death.  
Other parameters are fixed based on values from the literature ([Parameters](./Parameters)).
The model is fitted to the data using sequential Monte Carlo, a particle filter ([Inference](./Inference)).

### Model outputs
We simulate the dynamics of SARS-CoV-2 in Illinois using the best-fit parameters from the model inference while incorporating uncertainty.
Simulations involving different public health interventions will be later uploaded to the [Forecasting](./Forecasting) directory.
The forecasts incorporate several types of uncertainty that contribute to variation:
* Uncertainty in precisely how many individuals will be infected, recover, or die each day (demographic stochasticity)
* Uncertainty in the final inferred transmission rate for each region.

Incorporating this uncertainty produces a range of potential epidemic trajectories.

### Caveats
The model's predictions will shift as new data on COVID-19 emerges from Illinois and around the world.
We will update the model to incorporate better assumptions about the underlying biology and epidemiology.
We will also continue to try different modeling approaches with different assumptions to see how they affect conclusions.
And we will continue to engage with other scientists and modelers, and learn from what they are doing (e.g., via the [MIDAS Network](https://midasnetwork.us/) and the [Mobility Data Network](https://www.covid19mobility.org/)). 
[This recent perspective](https://www.nejm.org/doi/full/10.1056/NEJMp2016822) highlights important issues to consider when evaluating mathematical models.

Some important assumptions of our model include:
* Infected individuals who become hospitalized no longer contribute to population transmission. Obviously, if insufficient PPE is available or a hospitalized case is not diagnosed in time, this could be a bad assumption.
* The geographic regions of Illinois are independent (i.e., meaningful transmission does not occur between regions).
* Changes in season have no effect on susceptibility or transmission. Although there probably is some seasonal variation in these factors, which leads to wintertime colds and flus in temperate populations, the magnitude of these effects has been a longstanding question in infectious disease biology. Pandemics in particular tend to violate typical patterns. We discuss the complexity of these seasonal factors in a [recent paper](https://science.sciencemag.org/content/early/2020/04/23/science.abb5659/tab-article-info).
* Many parameters in our model are fixed based on existing literature (see [Parameters](./Parameters)).

### More information

Full scientific methods (technical documentation) are in preparation.
Pull requests can be submitted directly.
Other correspondence should be sent to cobey@uchicago.edu.

### References
1. King AA, Nguyen D and Ionides EL (2015) Statistical inference for partially observed Markov processes via the R package pomp. arXiv preprint arXiv:1509.00503.
2. Holmdahl, I. and Buckee, C (2020) Wrong but Useful — What Covid-19 Epidemiologic Models Can and Cannot Tell Us. <i>New England Journal of Medicine.</i> 10.1056/NEJMp2016822
3. Cobey, S (2020) Modelling infectious disease dynamics. <i>Science</i> eabb5659


### License for forecast data
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />Data in Forecasting/forecast_hub_projections/ are licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

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
