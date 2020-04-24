# Parameter values

Here, we show the values of parameters we fixed based on current literature on covid-19.


# Sources

## Parameters derived from the serial interval and incubation period: &zeta;<sub>symp.</sub> and &sigma;

Best estimates place the mean incubation period of the virus at 6 days (Backer, Lauer, Tindale, Linton) and the mean serial interval at 4.25 days (Du, Tindale, Nishura). Assuming that cases become infectious 0.5-1 days prior to transmission, we can use these estimates to calculate the time between infection and infectiousness (&sigma;<sup>-1</sup>) and the time from infectiousness to symptom onset (&zeta;<sub>symp.</sub><sup>-1</sup>).

## Rate of recovery for mildly symptomatic infections: &gamma;<sub>m</sub>

Based on viral isolation data from patients with mild disease, we estiamte that infectiousness of mildly symptomatic people wanes after 7.5 days (Wolfel).

## Rate of recovery for asymptomatic infections: &eta;

We assumed that the duration of infectiousness is similar in mild and asymptomatic cases and therefore estimate that the total time to recovery for asymptomatic people is 10 days.

## Fraction of infections that are asymptomatic: &rho;

## Fraction of presymptomatic infections that become severe: &rho;

## Test sensitivity and fraction of hospitalized infections that are tested: &theta;<sub>test</sub> and &nu;<sub>hosp.</sub>


# Parameters affecting hospitalization

Parameters that determine how hospitalized people are partitioned (&psi;<sub>1,2,3,4</sub>) and the rates at which they move through the hospital to death or recovery (&zeta;<sub>crit.</sub>, &gamma;<sub>hosp.</sub>, &mu;<sub>crit.</sub>, &mu;<sub>hosp.</sub>)were based on current data from a hospital system in Illinois. 
