# Parameter values

Here, we show the values of parameters we fixed based on current literature on covid-19.


# Sources

## Parameters derived from the serial interval and incubation period: &zeta;<sub>symp.</sub> and &sigma;

Best estimates place the mean incubation period of the virus at 6 days<sup>1,2,3,4</sup> and the mean serial interval at 4.25 days<sup>3,5,6</sup>. Assuming that cases become infectious 0.5-1 days prior to transmission, we can use these estimates to calculate the time between infection and infectiousness (&sigma;<sup>-1</sup>) and the time from infectiousness to symptom onset (&zeta;<sub>symp.</sub><sup>-1</sup>).

## Rate of recovery for mildly symptomatic infections: &gamma;<sub>m</sub>

Based on viral isolation data from patients with mild disease, we estiamte that infectiousness of mildly symptomatic people wanes after 7.5 days (Wolfel).

## Rate of recovery for asymptomatic infections: &eta;

We assumed that the duration of infectiousness is similar in mild and asymptomatic cases and therefore estimate that the total time to recovery for asymptomatic people is 10 days.

## Fraction of infections that are asymptomatic: &rho;

## Fraction of presymptomatic infections that become severe: &rho;

## Test sensitivity and fraction of hospitalized infections that are tested: &theta;<sub>test</sub> and &nu;<sub>hosp.</sub>


# Parameters affecting hospitalization

Parameters that determine how hospitalized people are partitioned (&psi;<sub>1,2,3,4</sub>) and the rates at which they move through the hospital to death or recovery (&zeta;<sub>crit.</sub>, &gamma;<sub>hosp.</sub>, &mu;<sub>crit.</sub>, &mu;<sub>hosp.</sub>)were based on current data from a hospital system in Illinois. 

# References

1.  Backer Jantien A, Klinkenberg Don, Wallinga Jacco. Incubation period of 2019 novel coronavirus (2019-nCoV) infections among travellers from Wuhan, China, 20–28 January 2020. Euro Surveill. 2020;25(5):pii=2000062. https://doi.org/10.2807/1560-7917.ES.2020.25.5.2000062
2. Lauer SA, Grantz KH, Bi Q, et al. The Incubation Period of Coronavirus Disease 2019 (COVID-19) From Publicly Reported Confirmed Cases: Estimation and Application. Ann Intern Med. 2020; [Epub ahead of print 10 March 2020]. doi: https://doi.org/10.7326/M20-0504
3. Lauren Tindale, Michelle Coombe, Jessica E Stockdale, Emma Garlock, Wing Yin Venus Lau, Manu Saraswat, Yen-Hsiang Brian Lee, Louxin Zhang, Dongxuan Chen, Jacco Wallinga, Caroline Colijn
medRxiv 2020.03.03.20029983; doi: https://doi.org/10.1101/2020.03.03.20029983 
4. Linton, N.M.; Kobayashi, T.; Yang, Y.; Hayashi, K.; Akhmetzhanov, A.R.; Jung, S.-M.; Yuan, B.; Kinoshita, R.; Nishiura, H. Incubation Period and Other Epidemiological Characteristics of 2019 Novel Coronavirus Infections with Right Truncation: A Statistical Analysis of Publicly Available Case Data. J. Clin. Med. 2020, 9, 538. 
5. Du Z, Xu X, Wu Y, Wang L, Cowling BJ, Ancel Meyers L. Serial interval of COVID-19 among publicly reported confirmed cases. Emerg Infect Dis. 2020 Jun. https://doi.org/10.3201/eid2606.200357
6. Nishiura H, Linton NM, Akhmetzhanov AR. Serial interval of novel coronavirus (COVID-19) infections. Int. J. Infect. Dis. 2020. doi: 10.1016/j.ijid.2020.02.060

