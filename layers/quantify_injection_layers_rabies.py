# -*- coding: utf-8 -*-
"""
Spyder Editor

Perform a Spearman correlation on the regions and layers inside the 
injection polygon - target-defined expeirments matched with wild type.
"""

import os
import pandas as pd
import json
from anatomy.anatomy_api import AnatomyApi
aapi = AnatomyApi()
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
ai_map = {value:key for key,value in ia_map.items()}
import platform
if platform.system() == 'Windows':
    path = r'2019 DMN\_new_figures\Figure_5'
else:
    path = r'2019 DMN/_new_figures/Figure_5'
dat = pd.read_csv(os.path.join(path, 'rabies_matches_for_td.csv'))
ss = aapi.get_summary_structure_data('id')
isocortex = st.get_structures_by_acronym(['Isocortex'])[0]

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
for isid in dat['id'].unique():
    relative_projs = []
    try:
        unionize = mcc.get_structure_unionizes([isid], is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    except:
        if platform.system() == 'Windows':
            unionize_path = os.path.join(r'2019 DMN\data_files\alternative_unionizes',
                                         'experiment_{0}'.format(str(isid)),
                                         'output.json')
        else:
            unionize_path = os.path.join(r'2019 DMN/data_files/alternative_unionizes',
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
        dat.loc[dat['id'] == isid, keys[ix]] = layer_projs/unionize[
            unionize['structure_id'] == isocortex['id']]['normalized_projection_volume'].sum()
    l23_projs = unionize[unionize['structure_id'].isin(L2_3strs)]['normalized_projection_volume'].sum()
    l5_projs = unionize[unionize['structure_id'].isin(L5strs)]['normalized_projection_volume'].sum()
    dat.loc[dat['id'] == isid, 'L23_5_ratio'] = l23_projs/l5_projs
dat.to_csv(os.path.join(path, 'rabies_injection_layers.csv'),
           index = False)
