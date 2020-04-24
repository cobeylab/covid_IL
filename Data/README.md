# Data for Illinois COVID-19 model

## Demography
`IL_population_EMS_regions_1-2.csv`, `IL_population_EMS_regions_3-6.csv`, and `IL_population_EMS_regions_7-11.csv` are the July 2018 Census population estimates in 10-year age bins for EMS regions 1-2, 3-6, and 7-11 ([IDPH](https://www.dph.illinois.gov/topics-services/emergency-preparedness-response/ems/preHospData)). We divided the state into these regions based on correlations in their observed incident death dynamics.

![alt text](EMS%20regions%20map.png)

## Contact data
To parameterize contact rates, we used age-structured contact matrices for the United States from Prem, et al. 2017<sup>1</sup>. There are four separate matrices describing contacts at home, school, work, and other locations. We reaggregated these matrices into 10-year age bins and provide them in `formatted_contacts_IL.RData`.

## Underreporting

To account for underreporting in deaths, we modified the approach of Weinberger et al.<sup>2</sup> to calculate the number of excess pneumonia and influenza (P & I) deaths observed in that timeframe in Illinois. `underreporting_PnI_IL.csv` contains estimates of the number of excess P & I deaths by week. By comparing these to the number of observed COVID-19 deaths, we were able to estimate the fraction of COVID-19 deaths that were not reported. Baed on this method, we find no underreporting of deaths after March 28. The final underreporting fractions we used in our model are in `frac_underreported.csv`. We also provide `nu_scaling.csv` which can scale reporting rates directly if desired. Currently, all scaling values are set to 1 (i.e., no additional scaling).

## References
1. Prem K, Cook AR, Jit M (2017) Projecting social contact matrices in 152 countries using contact surveys and demographic data. PLOS Computational Biology 13(9): e1005697. https://doi.org/10.1371/journal.pcbi.1005697
2. Daniel Weinberger, Ted Cohen, Forrest Crawford, Farzad Mostashari, Don Olson, Virginia E Pitzer, Nicholas G Reich, Marcus Russi, Lone Simonsen, Annie Watkins, Cecile Viboud. Estimating the early death toll of COVID-19 in the United States. medRxiv 2020.04.15.20066431; doi: https://doi.org/10.1101/2020.04.15.20066431 
