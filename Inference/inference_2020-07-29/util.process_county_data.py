import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import pandas as pd

def change_county_name(cty):
    if cty in ['Chicago', 'Cook', 'Suburban Cook']:
        new_cty = 'COOK'
    else:
        new_cty = cty.upper()
    return new_cty 

regions = pd.read_csv('/project2/cobey/covid_IL/Data/counties_census.csv')
regions['county'] = [cty.upper() for cty in regions.county]
regions.index = regions.county

idph = pd.read_csv('/project2/cobey/covid-civis//raw_data/case_data_public/idph_public_by_county.csv')
chicago = idph[idph.county =='Chicago'].copy()
chicago['region'] = 'chicago'
chicago['test_date'] = pd.to_datetime(chicago.test_date)
chicago = chicago.groupby(['region', 'test_date']).sum().reset_index()
chicago['inc_deaths'] = chicago.groupby('region').diff()['deaths']

idph['county'] = idph.county.apply(change_county_name)
idph = idph[~idph.county.isin(['UNASSIGNED', 'OUT OF STATE', 'ILLINOIS'])]
idph['region'] = idph.county.apply(lambda x: regions.loc[x, 'idph_region'])
idph['test_date'] = pd.to_datetime(idph.test_date)

incident = idph.groupby(['region', 'test_date']).sum().reset_index()
incident['inc_deaths'] = incident.groupby('region').diff()['deaths']

incident = incident.append(chicago)

for region, rdf in incident.groupby('region'):
    rdf.plot(x='test_date', y='inc_deaths')
    plt.title(region)

incident.dropna().to_csv('./incident_idph_regions.csv', index=False)
