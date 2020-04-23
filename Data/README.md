# Demography
`IL_population_EMS_regions_1-2.csv`, `IL_population_EMS_regions_3-6.csv`, and `IL_population_EMS_regions_7-11.csv` are the July 2018 Census population estimates in 10-year age bins for EMS regions 1-2, 3-6, and 7-11 ([IDPH](https://www.dph.illinois.gov/topics-services/emergency-preparedness-response/ems/preHospData)). We divided the state into these regions based on correlations in their observed incident death dynamics.

![alt text](EMS%20regions%20map.png)


# Contact data
To parameterize contact rates, we used age-structured contact matrices for the United States from [Prem, et al. 2017. PLoS Comp. Bio.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697). There are four separate matrices describing contacts at home, school, work, and other locations. We reaggregated these matrices into 10-year age bins and provide them in `formatted_contacts_IL.RData`.
