# -*- coding: utf-8 -*-
"""
Spyder Editor

Perform a Spearman correlation on the regions and layers in the injection site poluygon
Wild type matched data
"""

import os
import pandas as pd
import scipy.stats as stats
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
ai_map = {value:key for key,value in ia_map.items()}

import platform
if platform.system() == 'Windows':
    path = r'c:\\Users\jenniferwh\Dropbox\DMN data\correlations\_final'
elif platform.system() == 'Darwin':
    path = r'/Users/jenniferwh/Dropbox/DMN data/correlations/_final'
dat = pd.read_csv(os.path.join(path, 'match_correlations_by_source_NPV_all_thresholded_1_5.csv'))

structs = st.get_structures_by_set_id(
        [667481440, 667481441, 667481445, 667481446, 667481449])
acronyms = [structure['acronym'] for structure in structs]
ids = [structure['id'] for structure in structs]
ids += [1121, 526, 20, 52, 543, 664, 92, 712, 139, 727, 28, 60, 743] #ENTl, ENTm
L1strs = [ia_map[structure] for structure in acronyms if '1' in structure]
L1strs += [1121, 526]
L2_3strs = [ia_map[structure] for structure in acronyms if '2/3' in structure]
L2_3strs += [20, 543, 52, 664] #ENT
L4strs = [ia_map[structure] for structure in acronyms if '4' in structure]
L4strs += [92, 712]
L5strs = [ia_map[structure] for structure in acronyms if '5' in structure]
L5strs += [139, 727]
L6strs = [ia_map[structure] for structure in acronyms if '6' in structure]
L6strs = [structure for structure in L6strs if '6b' not in ai_map[structure]] #6b is unreliable
L6strs += [28, 743]

# remove L6b structures from quantified injection area. There should not be any but just in case.
L6bstrs = [structure for structure in L6strs if '6b' in ai_map[structure]]
structs = [structure for structure in structs if structure not in L6bstrs]

for isid in dat['match_A'].unique():
    anchor_unionize = mcc.get_structure_unionizes([isid], is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    inj_volume = anchor_unionize['volume'].sum()
    anchor_unionize['relative_projection_volume'] = anchor_unionize[
            'projection_volume'] * anchor_unionize['volume']/ inj_volume
    anchor_unionize = anchor_unionize[['experiment_id',
                                       'relative_projection_volume',
                                       'structure_id']]
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L1strs), 'structure_id'] = 1
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L2_3strs), 'structure_id'] = 2
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L4strs), 'structure_id'] = 4
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L5strs), 'structure_id'] = 5
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L6strs), 'structure_id'] = 6
    anchor_unionize = anchor_unionize.groupby(['experiment_id', 'structure_id']).sum().reset_index()
    
    matches = dat[dat['match_A'] == isid]['match_B'].unique()
    test_unionizes = mcc.get_structure_unionizes(matches, is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    test_inj_volumes = test_unionizes.groupby('experiment_id')['volume'].sum()
    normalized_test_unionizes = pd.DataFrame(columns = test_unionizes.columns)
    normalized_test_unionizes['relative_projection_volume'] = np.nan
    for match_id in test_unionizes['experiment_id'].unique():
        match_unionize = test_unionizes[test_unionizes['experiment_id'] == match_id]
        match_unionize['relative_projection_volume'] = match_unionize[
                'projection_volume'] * match_unionize['volume'] / test_inj_volumes[match_id]
        normalized_test_unionizes = pd.concat([normalized_test_unionizes, match_unionize])
    test_unionizes = normalized_test_unionizes[['experiment_id',
                                     'relative_projection_volume',
                                     'structure_id',]]
    test_unionizes.loc[test_unionizes['structure_id'].isin(L1strs), 'structure_id'] = 1
    test_unionizes.loc[test_unionizes['structure_id'].isin(L2_3strs), 'structure_id'] = 2
    test_unionizes.loc[test_unionizes['structure_id'].isin(L4strs), 'structure_id'] = 4
    test_unionizes.loc[test_unionizes['structure_id'].isin(L5strs), 'structure_id'] = 5
    test_unionizes.loc[test_unionizes['structure_id'].isin(L6strs), 'structure_id'] = 6
    test_unionizes = test_unionizes.groupby(['experiment_id', 'structure_id']).sum().reset_index()
    
    for isidB in test_unionizes['experiment_id'].unique():
        unionizes = pd.concat([anchor_unionize, 
                              test_unionizes[
                                      test_unionizes['experiment_id'] == isidB]], 
        join = 'outer')
        for structure in unionizes['structure_id'].unique():
            if len(unionizes[
                    (unionizes['experiment_id'] == isid) & 
                    (unionizes['structure_id'] == structure)]) == 0:
                unionizes = unionizes.append({'experiment_id': isid, 
                                  'relative_projection_volume': np.nan,
                                  'structure_id': structure}, ignore_index = True)
            if len(unionizes[(unionizes['experiment_id'] == isidB) &
                             (unionizes['structure_id'] == structure)]) == 0:
                unionizes = unionizes.append({'experiment_id': isidB, 
                                  'relative_projection_volume': np.nan,
                                  'structure_id': structure}, ignore_index = True)
        unionizes.sort_values(by=['structure_id'], inplace = True)
        x, y = unionizes[unionizes['experiment_id'] == isid]['relative_projection_volume'].values, unionizes[
                unionizes['experiment_id'] == isidB]['relative_projection_volume'].values
        nas = np.logical_or(np.isnan(x), np.isnan(y))
        try:
            pearsonr, _ = stats.pearsonr(x[~nas], y[~nas])
        except:
            print('pearsonr failure', isid, isidB)
        dat.loc[(dat['match_A'] == isid) & (dat['match_B'] == isidB),
                'inj_corr'] = pearsonr

# Change layers slightly for entorhinal cortex
L2strs = [ia_map[structure] for structure in acronyms if '2/3' in structure]
L2strs += [20, 543]
L3strs = [52, 664]

ENTdat = dat[dat['source'].isin(['ENTm', 'ENTl'])]
for isid in ENTdat['match_A'].unique():
    anchor_unionize = mcc.get_structure_unionizes([isid], is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    inj_volume = anchor_unionize['volume'].sum()
    anchor_unionize['relative_projection_volume'] = anchor_unionize[
            'projection_volume'] * anchor_unionize['volume']/ inj_volume
    anchor_unionize = anchor_unionize[['experiment_id',
                                       'relative_projection_volume',
                                       'structure_id']]
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L1strs), 'structure_id'] = 1
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L2strs), 'structure_id'] = 2
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L3strs), 'structure_id'] = 3
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L5strs), 'structure_id'] = 5
    anchor_unionize.loc[anchor_unionize['structure_id'].isin(L6strs), 'structure_id'] = 6
    anchor_unionize = anchor_unionize.groupby(['experiment_id', 'structure_id']).sum().reset_index()
    
    matches = dat[dat['match_A'] == isid]['match_B'].unique()
    test_unionizes = mcc.get_structure_unionizes(matches, is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    test_inj_volumes = test_unionizes.groupby('experiment_id')['volume'].sum()
    normalized_test_unionizes = pd.DataFrame(columns = test_unionizes.columns)
    normalized_test_unionizes['relative_projection_volume'] = np.nan
    for match_id in test_unionizes['experiment_id'].unique():
        match_unionize = test_unionizes[test_unionizes['experiment_id'] == match_id]
        match_unionize['relative_projection_volume'] = match_unionize[
                'projection_volume'] * match_unionize['volume'] / test_inj_volumes[
                        match_id]
        normalized_test_unionizes = pd.concat([normalized_test_unionizes, match_unionize])
    test_unionizes = normalized_test_unionizes[['experiment_id',
                                     'relative_projection_volume',
                                     'structure_id',]]
    test_unionizes.loc[test_unionizes['structure_id'].isin(L1strs), 'structure_id'] = 1
    test_unionizes.loc[test_unionizes['structure_id'].isin(L2strs), 'structure_id'] = 2
    test_unionizes.loc[test_unionizes['structure_id'].isin(L3strs), 'structure_id'] = 3
    test_unionizes.loc[test_unionizes['structure_id'].isin(L5strs), 'structure_id'] = 5
    test_unionizes.loc[test_unionizes['structure_id'].isin(L6strs), 'structure_id'] = 6
    test_unionizes = test_unionizes.groupby(['experiment_id', 'structure_id']).sum().reset_index()
    
    for isidB in test_unionizes['experiment_id'].unique():
        unionizes = pd.concat([anchor_unionize, 
                              test_unionizes[
                                      test_unionizes['experiment_id'] == isidB]], 
        join = 'outer')
        for structure in unionizes['structure_id'].unique():
            if len(unionizes[
                    (unionizes['experiment_id'] == isid) & 
                    (unionizes['structure_id'] == structure)]) == 0:
                unionizes = unionizes.append({'experiment_id': isid, 
                                  'relative_projection_volume': np.nan,
                                  'structure_id': structure}, ignore_index = True)
            if len(unionizes[(unionizes['experiment_id'] == isidB) &
                             (unionizes['structure_id'] == structure)]) == 0:
                unionizes = unionizes.append({'experiment_id': isidB, 
                                  'relative_projection_volume': np.nan,
                                  'structure_id': structure}, ignore_index = True)
        unionizes.sort_values(by=['structure_id'], inplace = True)
        x, y = unionizes[unionizes['experiment_id'] == isid]['relative_projection_volume'].values, unionizes[
                unionizes['experiment_id'] == isidB]['relative_projection_volume'].values
        nas = np.logical_or(np.isnan(x), np.isnan(y))
        pearsonr, _ = stats.pearsonr(x[~nas], y[~nas])
        dat.loc[(dat['match_A'] == isid) & (dat['match_B'] == isidB),
                'inj_corr'] = pearsonr
     
dat.to_csv(os.path.join(path, 'match_correlations_by_source_NPV_all_thresh_inj_corr.csv'), 
           index  = False)