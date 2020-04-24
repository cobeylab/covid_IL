# Forecasting COVID-19 dynamics for the state of Illinois

*Contributors (alphabetically listed): Phil Arevalo, Ed Baskerville, Spencer Carran, Sarah Cobey, Katelyn Gostic, Lauren McGough, Sylvia Ranjeva, & Frank Wen 

### Introduction 
We developed a model to infer key aspects of COVID-19 transmission in Illinois, and to forecast community spread, hospital and ICU burden, and mortality under current and hypothetical public health interventions. 

### Model overview
We introduce an adaptation of an age-structured SEIR model (Figure 1) that accounts for asymptomatic infection and mortality. Individuals are infected, but not infectious upon entry into the exposed class (E). Individuals that enter the asymptomatic class (A) can infect others in the community but will eventually clear the infection without ever showing symptoms of disease, while those that enter the presymptomatic infectious class (P) will progress to symptomatic disease. Symptomatic infections are divided into mild cases, I<sub>M</sub>, which will resolve without hospital attention, and severe cases, I<sub>S</sub> that will require hospitalization. The class I<sub>H</sub> refers to severe cases that are hospitalized but not in the ICU, and the I<sub>C</sub> class refers to severe cases in the ICU. We explicitly track deaths from individuals in the hospitalized and ICU classes (I<sub>H4</sub> and I<sub>C</sub>, respectively). We account for deaths that occur outside of the hospital in our simulations via a scaling factor informed by the model parameters and by data on outside-of-hospital deaths from New York City (see [Parameters](./Parameters)).

![Figure 1](model_diagram.png)
Rate parameters determine the amount of time the average infected person spends in a given compartment. Probabilities determine what fraction of infected people follow specific transition paths between compartments.
To incorporate demographic stochasticity, the model is implemented in a sub-day discrete-time approximation of a continuous-time stochastic framework.

### Data
We model the state-level dynamics of COVID-19 via three distinct geographic regions, as described in [Data](./Data). The model is fitted to in-hospital deaths reported by the New York Times from March 15 to March 24, 2020 and to in-hospital deaths reported by the Illinois Department of Publc Health from March 24, 2020 onwards. 

### Inference
For each region in Illinois, we infer the transmission rate of COVID-19 before and after public health interventions, as well as the number of individuals infected at the beginning of the simulations (on March 1, 2020). We fix all other parameters based on values chosen from the literature, as described in the [Parameters](./Parameters) directory. 
We fitted the model to the data using maximum likelihood methods for partially observed Markov process (POMP) models (details in [Inference](./Inference)).

### Model outputs
We simulate the dynamics of COVID-19 in Illinois using the best-fit parameters from the model inference.
Our model can forecast the effects of public health interventions of varying forms, strengths, and durations. 
These model forecasts will be provided in the [Forecasting](./Forecasting) directory. 

### Public health interventions 
We incorporate interventions as scaling factors on the transmission rate for all infected individuals, and as reductions in  age- and setting-specific contact rates. Our baseline model scenario reflects the enactment of "shelter in place" interventions in Illinois beginning on March 16, 2020, extended indefinitely ("indefinite"). The shelter in place scenario involves an inferred reduction in transmission rate, a 100% reduction in at-school contacts, a 40% reduction in at-work contacts, and a 50% reduction in all other contacts occurring outside of the home. We also consider two hypothetical scenarios:

1. Lifting shelter in place on May 1, 2020 ("lifted")
2. The absence of shelter in place ("never")
  
### References
1. King AA, Nguyen D and Ionides EL (2015) Statistical inference for partially observed Markov processes via the R package pomp. arXiv preprint arXiv:1509.00503.

-------------------------------------------

Distributed under the xx --- insert CC license---.
