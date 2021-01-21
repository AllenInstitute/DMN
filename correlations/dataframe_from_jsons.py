# -*- coding: utf-8 -*-
"""
Created on Fri Sep 07 17:48:30 2018

@author: JENNIFERWH
"""
# ",".join(map( lambda x: str(x), lst2))

import os
import json
import pandas as pd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi

aapi = AnatomyApi() 
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
ss_ids = aapi.get_summary_structure_data('id')
ss_acronyms = aapi.get_summary_structure_data('acronym')
st = mcc.get_structure_tree()
ai_map = st.get_id_acronym_map()
ia_map = {value: key for key, value in ai_map.items()}

#Compile wt data
path = r'cluster_code\correlations\output\NPV'
df = pd.DataFrame()
for filename in os.listdir(path):
    if filename != 'distance_corrected' and filename != 'good backup':
        jsondat = os.path.join(path, filename)
        with open(jsondat, 'r') as data_file:    
            data = json.load(data_file)
            dat = pd.DataFrame(data)
            df = pd.concat([df, dat])

for experiment in df['image_series_id']:
    unionize = mcc.get_experiment_structure_unionizes(experiment_id = experiment)
    if len(unionize) > 0:
        unionize = unionize[(unionize['is_injection'] == True) &
                            (unionize['hemisphere_id'] == 3)]
    else:
        unionize_path = os.path.join(r'DMN_paper\alternative_unionizes',
                                         'experiment_{0}'.format(str(experiment)), 
                                         'output.json') #new data not online yet
        with open(unionize_path, 'r') as jsonfile:
            unionize = json.load(jsonfile)
        unionize = pd.DataFrame(unionize)
        unionize = unionize[unionize['threshold'] == 0]
        unionize.drop_duplicates(inplace = True)
        unionize.rename(columns = {'projection_volume_above_threshold': 'projection_volume',
                                   'normalized_projection_volume_above_threshold': 
                                       'normalized_projection_volume',
                                       'image_series_id': 'experiment_id'}, 
            inplace = True)
        unionize = unionize[unionize['hemisphere_id'] == 3]
    inj_size = unionize[unionize['structure_id'] == 997]['projection_volume'].values[0]
    ss_dat = unionize[unionize['structure_id'].isin(ss_ids)]
    primary = df[df['image_series_id'] == experiment]['source'].values[0]
    percent_primary = ss_dat[ss_dat['structure_id'] == ai_map[primary]]['projection_volume'].values[0]/inj_size
    if len(ss_dat) > 1:
        secondary = ss_dat.sort_values(
            by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[1]
        percent_secondary = ss_dat.sort_values(
            by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[1]/inj_size
    df.loc[df['image_series_id'] == experiment, 'td_secondary_source'] = ia_map[secondary]
    df.loc[df['image_series_id'] == experiment, 'td_percent_secondary'] = percent_secondary

    df.loc[df['image_series_id'] == experiment, 'td_injection_size'] = inj_size
    df.loc[df['image_series_id'] == experiment, 'td_primary_source'] = primary
    df.loc[df['image_series_id'] == experiment, 'td_percent_primary'] = percent_primary

alldat = pd.DataFrame(mcc.get_experiments())
for experiment in df['match_id'].unique():
    unionize = mcc.get_experiment_structure_unionizes(experiment_id = experiment)
    unionize = unionize[(unionize['is_injection'] == True) &
                        (unionize['hemisphere_id'] == 3)]
    inj_size = unionize[unionize['structure_id'] == 997]['projection_volume'].values[0]
    ss_dat = unionize[unionize['structure_id'].isin(ss_ids)]
    primary = ss_dat.sort_values(
        by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[0]
    percent_primary = ss_dat.sort_values(
        by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[0]/inj_size
    try:
        secondary = ss_dat.sort_values(
            by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[1]
        percent_secondary = ss_dat.sort_values(
            by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[1]/inj_size
        df.loc[df['match_id'] == experiment, 'match_secondary_source'] = ia_map[secondary]
        df.loc[df['match_id'] == experiment, 'match_percent_secondary'] = percent_secondary
    except:
        df.loc[df['match_id'] == experiment, 'match_secondary_source'] = np.nan
        df.loc[df['match_id'] == experiment, 'match_percent_secondary'] = np.nan
    df.loc[df['match_id'] == experiment, 'match_injection_size'] = inj_size
    df.loc[df['match_id'] == experiment, 'match_primary_source'] = ia_map[primary]
    df.loc[df['match_id'] == experiment, 'match_percent_primary'] = percent_primary
    df.loc[df['match_id'] == experiment, 'transgenic_line'] = alldat[alldat['id'] == experiment]['transgenic_line'].values[0]
    if alldat[alldat['id'] == experiment]['product_id'].values == 36:
        df.loc[df['match_id'] == experiment, 'Virus'] = 'SypEGFP'
    else:
        df.loc[df['match_id'] == experiment, 'Virus'] = 'EGFP'

df['same_primary'] = df['td_primary_source'] == df['match_primary_source']
df['same_secondary'] = df['td_secondary_source'] == df['match_secondary_source']
df['same secondary for <60% primary'] = (df['td_percent_primary'] < 0.6)  | (
        df['match_percent_primary'] < 0.6) & (df['same_primary']) & (df['same_secondary'])
df.loc[df['same_primary'] == True, 'same secondary for <60% primary'] = np.nan
df['injection_size_ratio'] = df['match_injection_size']/df['td_injection_size']   
df.to_csv(r'DMN data\correlations\td_wt_cre_matched_correlations_NPV_dice.csv', index = False)
