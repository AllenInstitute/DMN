# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import os
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi

import platform
if platform.system() == 'Darwin':
    path = '/Users/jenniferwh/Dropbox/DMN data/correlations'
elif platform.system() == 'Windows':
    path = r'C:\\Users\jenniferwh\Dropbox\DMN data\correlations'
    
aapi = AnatomyApi()
data = pd.read_csv(os.path.join(path, 'cla_correlations.csv'))

mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
ss_ids = aapi.get_summary_structure_data('id')
ss_acronyms = aapi.get_summary_structure_data('acronym')
st = mcc.get_structure_tree()
ai_map = st.get_id_acronym_map()
ia_map = {value: key for key, value in ai_map.items()}

for experiment in data['image_series_id']:
    unionize = mcc.get_experiment_structure_unionizes(experiment_id = experiment, 
                                           hemisphere_ids = [3],
                                           is_injection = True)
    inj_size = unionize[unionize['structure_id'] == 997]['projection_volume'].values[0]
    ss_dat = unionize[unionize['structure_id'].isin(ss_ids)]
    primary = ss_dat.sort_values(by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[0]
    percent_primary = ss_dat.sort_values(
        by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[0]/inj_size
    if len(ss_dat) > 1:
        secondary = ss_dat.sort_values(
            by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[1]
        percent_secondary = ss_dat.sort_values(
            by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[1]/inj_size
        data.loc[data['image_series_id'] == experiment, 'td_secondary_source'] = ia_map[secondary]
        data.loc[data['image_series_id'] == experiment, 'td_percent_secondary'] = percent_secondary
    data.loc[data['image_series_id'] == experiment, 'td_percent_primary'] = percent_primary
    
for experiment in data['match_id'].unique(): 
    unionize = mcc.get_experiment_structure_unionizes(experiment_id = experiment, 
                                           hemisphere_ids = [3],
                                           is_injection = True)
    inj_size = unionize[unionize['structure_id'] == 997]['projection_volume'].values[0]
    ss_dat = unionize[unionize['structure_id'].isin(ss_ids)]
    primary = ss_dat.sort_values(by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[0]
    percent_primary = ss_dat.sort_values(
        by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[0]/inj_size
    if len(ss_dat) > 1:
        secondary = ss_dat.sort_values(
            by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[1]
        percent_secondary = ss_dat.sort_values(
            by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[1]/inj_size
        
        data.loc[data['match_id'] == experiment, 'match_secondary_source'] = ia_map[secondary]
        data.loc[data['match_id'] == experiment, 'match_percent_secondary'] = percent_secondary
    data.loc[data['match_id'] == experiment, 'match_primary_source'] = ia_map[primary]
    data.loc[data['match_id'] == experiment, 'match_percent_primary'] = percent_primary
    print(experiment)

data.to_csv(os.path.join(path, 'cla_correlations.csv'), index = False)