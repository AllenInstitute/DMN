# -*- coding: utf-8 -*-
"""
Spyder Editor

Perform a Spearman correlation on the regions and layers inside the 
injection polygon - target-defined expeirments matched with wild type.
"""

import os
import pandas as pd
import scipy.stats as stats
import numpy as np
import json
import nrrd
from anatomy.anatomy_api import AnatomyApi
aapi = AnatomyApi()
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
ai_map = {value:key for key,value in ia_map.items()}

path = r'/Users/jenniferwh/Dropbox/DMN data/correlations/_final'
dat = pd.read_csv(os.path.join(path, 'good_td_wt_correlations.csv'))

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

# remove L6b structures from quantified injection area
L6bstrs = [structure for structure in L6strs if '6b' in ai_map[structure]]
structs = [structure for structure in structs if structure not in L6bstrs]

for isid in dat['image_series_id'].unique():
    try:
        anchor_unionize = mcc.get_structure_unionizes([isid], is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    except:
        unionize_path = os.path.join(r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/data_files/alternative_unionizes',
                                 'experiment_{0}'.format(str(isid)), 
                                 'output.json') #new data not online yet
        if not os.path.isfile(unionize_path): #subcortical, non-hippocampal injections
            anchor_unionize = pd.DataFrame({
                    'experiment_id': isid,
                    'volume': 0,
                    'projection_volume': 0,
                    'structure_id': ids,
                    'hemisphere_id': 3})
        else:
            with open(unionize_path, 'r') as jsonfile:
                unionize_dat = json.load(jsonfile)
            unionize_dat = pd.DataFrame(unionize_dat)
            unionize_dat = unionize_dat[unionize_dat['threshold'] == 0]
            unionize_dat.drop_duplicates(inplace = True)
            unionize_dat.rename(columns = {'projection_volume_above_threshold': 'projection_volume',
                                       'normalized_projection_volume_above_threshold': 
                                           'normalized_projection_volume',
                                           'image_series_id': 'experiment_id'}, 
                inplace = True)
            anchor_unionize = unionize_dat[(unionize_dat['is_injection'] == True) &
                                           (unionize_dat['structure_id'].isin(ids)) &
                                           (unionize_dat['hemisphere_id'] == 3)]
            storage_directory = os.path.join(aapi.get_storage_directory(isid), 'grid')
            data_maskA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                               'data_mask_25.nrrd'))
            injA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                             'injection_density_25.nrrd'))
            injA = np.multiply(injA, data_maskA)
            injA_mask = injA[np.where(injA > 0.001)] = 1
            for structure in anchor_unionize['structure_id'].unique():
                structure_mask = mcc.get_structure_mask(structure)
                anchor_unionize.loc[anchor_unionize['structure_id'] == structure,
                                      'volume'] = np.sum(injA_mask[np.where(structure_mask)])
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
    
    matches = dat[dat['image_series_id'] == isid]['match_id'].unique()
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
        pearsonr, _ = stats.pearsonr(x[~nas], y[~nas])
        dat.loc[(dat['image_series_id'] == isid) & (dat['match_id'] == isidB),
                'inj_corr'] = pearsonr

# Change layers slightly for entorhinal cortex
L2strs = [ia_map[structure] for structure in acronyms if '2/3' in structure]
L2strs += [20, 543]
L3strs = [52, 664]

#%%
ENTdat = dat[dat['source'].isin(['ENTm', 'ENTl'])]
print(len(ENTdat))
for isid in ENTdat['image_series_id'].unique():
    try:
        anchor_unionize = mcc.get_structure_unionizes([isid], is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    except:
        unionize_path = os.path.join(r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/data_files/alternative_unionizes',
                                     'experiment_{0}'.format(str(isid)),
                                     'output.json') #new data not online yet
    else:
        with open(unionize_path, 'r') as jsonfile:
            unionize_dat = json.load(jsonfile)
        unionize_dat = pd.DataFrame(unionize_dat)
        unionize_dat = unionize_dat[unionize_dat['threshold'] == 0]
        unionize_dat.drop_duplicates(inplace = True)
        unionize_dat.rename(columns = {'projection_volume_above_threshold': 'projection_volume',
                                   'normalized_projection_volume_above_threshold': 
                                       'normalized_projection_volume',
                                       'image_series_id': 'experiment_id'}, 
            inplace = True)
        anchor_unionize = unionize_dat[(unionize_dat['is_injection'] == True) &
                                       (unionize_dat['structure_id'] == ids) &
                                       (unionize_dat['hemisphere_id'] == 3)]
        storage_directory = os.path.join(aapi.get_storage_path(isid), 'grid')
        data_maskA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                           'data_mask_25.nrrd'))
        injA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                         'injection_density_25.nrrd'))
        injA = np.multiply(injA, data_maskA)
        injA_mask = injA[np.where(injA > 0.001)] = 1
        for structure in anchor_unionize['structure_id'].unique():
            structure_mask = mcc.get_structure_mask(structure)
            anchor_unionize.loc[anchor_unionize['structure_id'] == structure,
                                  'volume'] = np.sum(injA_mask[np.where(structure_mask)])
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
    
    matches = dat[dat['image_series_id'] == isid]['match_id'].unique()
    test_unionizes = mcc.get_structure_unionizes(matches, is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    test_inj_volumes = test_unionizes.groupby('experiment_id')['volume'].sum()
    normalized_test_unionizes = pd.DataFrame(columns = test_unionizes.columns)
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
        dat.loc[(dat['image_series_id'] == isid) & (dat['match_id'] == isidB),
                'inj_corr'] = pearsonr

dat.to_csv(os.path.join(path, 'good_td_wt_correlations_with_inj_corr.csv'), 
           index  = False)