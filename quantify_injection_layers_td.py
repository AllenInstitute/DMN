# -*- coding: utf-8 -*-
"""
Spyder Editor

Perform a Spearman correlation on the regions and layers inside the 
injection polygon - target-defined expeirments matched with wild type.
"""

import os
import pandas as pd
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
import platform
if platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN'
else:
    path = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/'
dat = pd.read_csv(os.path.join(path, 'target_defined_dataset.csv'), engine = 'python')
dat = dat[dat['include'] == 'yes']
ENTdat = dat[dat['source'].isin(['ENTm', 'ENTl'])]
ss = aapi.get_summary_structure_data('id')
isocortex = st.get_structures_by_acronym(['Isocortex'])[0]
iso = st.descendant_ids([isocortex['id']])[0]
iso = [structure for structure in iso if structure in ss]
iso_strs = [ai_map[structure] for structure in iso]
dat = dat[dat['source'].isin(iso_strs)]

structs = st.get_structures_by_set_id(
        [667481440, 667481441, 667481445, 667481446, 667481449])
acronyms = [structure['acronym'] for structure in structs]
ids = [structure['id'] for structure in structs]
#ids += [1121, 526, 20, 52, 543, 664, 92, 712, 139, 727, 28, 60, 743] #ENTl, ENTm
L1strs = [ia_map[structure] for structure in acronyms if '1' in structure]
L2_3strs = [ia_map[structure] for structure in acronyms if '2/3' in structure]
L4strs = [ia_map[structure] for structure in acronyms if '4' in structure]
L5strs = [ia_map[structure] for structure in acronyms if '5' in structure]
L6strs = [ia_map[structure] for structure in acronyms if '6' in structure]
L6strs = [structure for structure in L6strs if '6b' not in ai_map[structure]] #6b is unreliable

# remove L6b structures from quantified injection area
L6bstrs = [structure for structure in L6strs if '6b' in ai_map[structure]]
ids = [structure for structure in ids if structure not in L6bstrs]
ids.append(isocortex['id'])
print([structure for structure in ids if structure not in L1strs+L2_3strs+L4strs+L5strs+L6strs])

# get relative volume of each layer
iso_mask, _ = mcc.get_structure_mask(isocortex['id'])
keys=['L1', 'L2/3', 'L4', 'L5', 'L6']

#%%
for isid in dat['image_series_id'].unique():
    relative_projs = []
    try:
        unionize = mcc.get_structure_unionizes([isid], is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    except:
        if platform.system() == 'Windows':
            unionize_path = os.path.join(r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\data_files\alternative_unionizes',
                                         'experiment_{0}'.format(str(isid)),
                                         'output.json')
        else:
            unionize_path = os.path.join(r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/data_files/alternative_unionizes',
                                 'experiment_{0}'.format(str(isid)), 
                                 'output.json') #new data not online yet
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
        unionize = unionize_dat[(unionize_dat['is_injection'] == True) &
                                       (unionize_dat['structure_id'].isin(ids)) &
                                       (unionize_dat['hemisphere_id'] == 3)]

    for ix, layer in enumerate([L1strs, L2_3strs, L4strs, L5strs, L6strs]):
        layer_projs = unionize[unionize['structure_id'].isin(layer)]['normalized_projection_volume'].sum()
        dat.loc[dat['image_series_id'] == isid, keys[ix]] = layer_projs/unionize[
            unionize['structure_id'] == isocortex['id']]['normalized_projection_volume'].sum()
    l23_projs = unionize[unionize['structure_id'].isin(L2_3strs)]['normalized_projection_volume'].sum()
    l5_projs = unionize[unionize['structure_id'].isin(L5strs)]['normalized_projection_volume'].sum()
    dat.loc[dat['image_series_id'] == isid, 'L23_5_ratio'] = l23_projs/l5_projs
dat.to_csv(os.path.join(path, '_new_figures', 'Figure_5', 'TD_injection_layers.csv'),
           index = False)

#%% This doesn't seem right. Not using for now
layer_volumes = dict()
for ix, layer in enumerate([L1strs, L2_3strs, L4strs, L5strs, L6strs]):
    big_mask = np.zeros(iso_mask.shape)
    for structure in layer:
        structure_mask, _ = mcc.get_structure_mask(structure)
        big_mask[np.where(structure_mask)] = 1
    layer_volumes[keys[ix]] = np.sum(big_mask)/np.sum(iso_mask)
    
#%%
# Change layers slightly for entorhinal cortex
L2strs = [ia_map[structure] for structure in acronyms if '2/3' in structure]
L2strs += [20, 543]
L3strs = [52, 664]

#%%
print(len(ENTdat))
L1strs += [1121, 526]
L2_3strs += [20, 543, 52, 664] #ENT
L4strs += [92, 712]
L5strs += [139, 727]
L6strs += [28, 743]
for isid in ENTdat['image_series_id'].unique():
    try:
        unionize2 = mcc.get_structure_unionizes([isid], is_injection = True,
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
        unionize2 = unionize_dat[(unionize_dat['is_injection'] == True) &
                                (unionize_dat['structure_id'].isin(ids)) &
                                (unionize_dat['hemisphere_id'] == 3)]
        storage_directory = os.path.join(aapi.get_storage_directory(isid))
        data_maskA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                           'data_mask_25.nrrd'))
        injA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                         'injection_density_25.nrrd'))
        injA = np.multiply(injA, data_maskA)
        injA_mask = injA_mask.copy()
        injA_mask[np.where(injA_mask > 0.001)] = 1
        for structure in unionize2['structure_id'].unique():
            structure_mask, _ = mcc.get_structure_mask(structure)
            unionize2.loc[unionize2['structure_id'] == structure,
                                  'volume'] = np.sum(injA_mask[np.where(structure_mask)])
    inj_volume = unionize2['volume'].sum()
    unionize2['relative_projection_volume'] = unionize2[
            'projection_volume'] * unionize2['volume']/ inj_volume
    unionize2 = unionize2[['experiment_id',
                         'relative_projection_volume',
                         'structure_id']]
    unionize2.loc[unionize2['structure_id'].isin(L1strs), 'structure_id'] = 1
    unionize2.loc[unionize2['structure_id'].isin(L2strs), 'structure_id'] = 2
    unionize2.loc[unionize2['structure_id'].isin(L3strs), 'structure_id'] = 3
    unionize2.loc[unionize2['structure_id'].isin(L5strs), 'structure_id'] = 5
    unionize2.loc[unionize2['structure_id'].isin(L6strs), 'structure_id'] = 6
    unionize2 = unionize2.groupby(['experiment_id', 'structure_id']).sum().reset_index()
    unionizes = pd.concat([unionize, unionize2])

unionizes.to_csv(os.path.join(path, '_new_figures', 'Figure_5', 'injection_layers_td_experiments.csv'), 
           index  = False)