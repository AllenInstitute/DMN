# -*- coding: utf-8 -*-
"""
Created on Fri Sep 07 17:48:30 2018

@author: JENNIFERWH
"""
# ",".join(map( lambda x: str(x), lst2))

import os
import json
import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi


#Compile wt data
path = r'cluster_code\correlations\output\matches_by_source\NPV'
df = pd.DataFrame()
for filename in os.listdir(path):
    print(filename)
    if filename not in ['all cre lines', 'with_non-matched_primaries',
                        'binarized', 'zeros_dropped', 'filtered',
                        'distance_corrected', 'good backup']:
        jsondat = os.path.join(path, filename)
        with open(jsondat, 'r') as data_file:    
            data = json.load(data_file)
            if len(data['pearson_correlation']) > 0:
                dat = pd.DataFrame(data)
                df = pd.concat([df, dat])

mcc = MouseConnectivityCache(manifest_file = 'connectivity\mouse_connectivity_manifest.json')
alldat = pd.DataFrame(mcc.get_experiments())
aapi = AnatomyApi()
ss_ids = aapi.get_summary_structure_data('id')
ss_acronyms = aapi.get_summary_structure_data('acronym')
st = mcc.get_structure_tree()
ai_map = st.get_id_acronym_map()
ia_map = {value: key for key, value in ai_map.items()}

for isid in df['match_A'].unique():
    unionize = mcc.get_experiment_structure_unionizes(experiment_id = isid, 
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
        df.loc[df['match_A'] == isid, 'match_A_secondary_source'] = ia_map[secondary]
        df.loc[df['match_B'] == isid, 'match_B_secondary_source'] = ia_map[secondary]
        df.loc[df['match_A'] == isid, 'match_A_percent_secondary'] = percent_secondary
        df.loc[df['match_B'] == isid, 'match_B_percent_secondary'] = percent_secondary
    df.loc[df['match_A'] == isid, 'match_A_injection_size'] = inj_size
    df.loc[df['match_B'] == isid, 'match_B_injection_size'] = inj_size
    df.loc[df['match_A'] == isid, 'match_A_primary_source'] = ia_map[primary]
    df.loc[df['match_B'] == isid, 'match_B_primary_source'] = ia_map[primary]
    df.loc[df['match_A'] == isid, 'match_A_percent_primary'] = percent_primary
    df.loc[df['match_B'] == isid, 'match_B_percent_primary'] = percent_primary
    df.loc[df['match_A'] == isid, 'transgenic_line_A'] = alldat[alldat['id'] == isid]['transgenic_line'].values[0]
    df.loc[df['match_B'] == isid, 'transgenic_line_B'] = alldat[alldat['id'] == isid]['transgenic_line'].values[0]
    
    if alldat[alldat['id'] == isid]['product_id'].values == 36:
        df.loc[df['match_A'] == isid, 'Virus_A'] = 'SypEGFP'
        df.loc[df['match_B'] == isid, 'Virus_B'] = 'SypEGFP'
    else:
        df.loc[df['match_A'] == isid, 'Virus_A'] = 'EGFP'
        df.loc[df['match_B'] == isid, 'Virus_B'] = 'EGFP'

for isid in df['match_B'].unique(): 
    if any(df[df['match_B'] == isid]['match_B_percent_primary'].isnull()):
        unionize = mcc.get_experiment_structure_unionizes(experiment_id = isid, 
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
            df.loc[df['match_B'] == isid, 'match_B_secondary_source'] = ia_map[secondary]
            df.loc[df['match_B'] == isid, 'match_B_percent_secondary'] = percent_secondary
        df.loc[df['match_B'] == isid, 'match_B_injection_size'] = inj_size
        df.loc[df['match_B'] == isid, 'match_B_primary_source'] = ia_map[primary]
        df.loc[df['match_B'] == isid, 'match_B_percent_primary'] = percent_primary
        df.loc[df['match_B'] == isid, 'transgenic_line_B'] = alldat[alldat['id'] == isid]['transgenic_line'].values[0]
        
        if alldat[alldat['id'] == isid]['product_id'].values == 36:
            df.loc[df['match_B'] == isid, 'Virus_B'] = 'SypEGFP'
        else:
            df.loc[df['match_B'] == isid, 'Virus_B'] = 'EGFP'

df2 = df.copy().reset_index()
# assign smaller injection to match A, larger to match B
for index, row in df2.iterrows():
    if(row['match_A_injection_size'] > row['match_B_injection_size']):
        df2.loc[index,['match_A','match_B']] = df2.loc[
                index,['match_B','match_A']].values
        df2.loc[index,['match_A_injection_size','match_B_injection_size']] = df2.loc[
                index,['match_B_injection_size','match_A_injection_size']].values
        df2.loc[index,['fraction of match covered by td injection','injection_overlap']] = df2.loc[
                index,['injection_overlap','fraction of match covered by td injection']].values
        df2.loc[index,['transgenic_line_A','transgenic_line_B']] = df2.loc[
                index,['transgenic_line_B','transgenic_line_A']].values
        df2.loc[index,['Virus_A','Virus_B']] = df2.loc[
                index,['Virus_B','Virus_A']].values
        df2.loc[index,['match_A_primary_source','match_B_primary_source']] = df2.loc[
                index,['match_B_primary_source','match_A_primary_source']].values
        df2.loc[index,['match_A_secondary_source','match_B_secondary_source']] = df2.loc[
                index,['match_B_secondary_source','match_A_secondary_source']].values
        df2.loc[index,['match_A_percent_primary','match_B_percent_primary']] = df2.loc[
                index,['match_B_percent_primary','match_A_percent_primary']].values
        df2.loc[index,['match_A_percent_secondary','match_B_percent_secondary']] = df2.loc[
                index,['match_B_percent_secondary','match_A_percent_secondary']].values
df2['same_primary'] = df2['match_A_primary_source'] == df2['match_B_primary_source']
df2['same_secondary'] = df2['match_A_secondary_source'] == df2['match_B_secondary_source']
df2['same secondary for <60% primary'] = (df2['match_A_percent_primary'] < 0.5)  | (
        df2['match_B_percent_primary'] < 0.5) & (df2['same_primary']) & (df2['same_secondary'])
                
df2.to_csv(r'DMN data\correlations\match_correlations_by_source_NPV_dice.csv', index = False)
