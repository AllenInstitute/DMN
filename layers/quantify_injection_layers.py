# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import pandas as pd
import scipy.stats as stats
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()

path = r'/Users/jenniferwh/Dropbox/DMN data/correlations'
dat = pd.read_csv(os.path.join(path, 'match_correlations_by_source_NPV_ipsi.csv'))

structs = st.get_structures_by_set_id(
        [667481440, 667481441, 667481445, 667481446, 667481449, 667481450])
acronyms = [structure['acronym'] for structure in structs]
ids = [structure['id'] for structure in structs]
ids += [1121, 526, 20, 52, 543, 664, 92, 712, 139, 727, 28, 60, 743] #ENTl, ENTm
L1strs = [ia_map[structure] for structure in acronyms if '1' in structure]
L1strs += [1121, 526]
L2_3strs = [ia_map[structure] for structure in acronyms if '2/3' in structure]
L2_3strs += [20, 52, 543, 664] #ENT
L4strs = [ia_map[structure] for structure in acronyms if '1' in structure]
L4strs += [92, 712]
L5strs = [ia_map[structure] for structure in acronyms if '5' in structure]
L5strs += [139, 727]
L6strs = [ia_map[structure] for structure in acronyms if '6' in structure]
L6strs += [28, 60, 743]

for isid in dat['match_A'].unique():
    anchor_unionize = mcc.get_structure_unionizes([isid], is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    anchor_unionize = anchor_unionize[['experiment_id',
                                       'projection_volume',
                                       'structure_id']]
    matches = dat[dat['match_A'] == isid]['match_B'].unique()
    test_unionizes = mcc.get_structure_unionizes(matches, is_injection = True,
                                       structure_ids = ids, hemisphere_ids = [3])
    test_unionizes = test_unionizes[['experiment_id',
                                     'projection_volume',
                                     'structure_id']]
    for isidB in matches:
        unionizes = pd.concat([anchor_unionize, 
                              test_unionizes[
                                      test_unionizes['experiment_id'] == isidB]], 
        join = 'outer')
        for structure in unionizes['structure_id'].unique():
            if len(unionizes[
                    (unionizes['experiment_id'] == isid) & 
                    (unionizes['structure_id'] == structure)]) == 0:
                unionizes = unionizes.append({'experiment_id': isid, 
                                  'projection_volume': 0,
                                  'structure_id': structure}, ignore_index = True)
            if len(unionizes[(unionizes['experiment_id'] == isidB) &
                             (unionizes['structure_id'] == structure)]) == 0:
                unionizes = unionizes.append({'experiment_id': isidB, 
                                  'projection_volume': 0,
                                  'structure_id': structure}, ignore_index = True)
        unionizes.sort_values(by=['structure_id'], inplace = True)
        spearmanr, _ = stats.spearmanr(unionizes[
                unionizes['experiment_id'] == isid]['projection_volume'], 
                    unionizes[unionizes['experiment_id'] == isidB]['projection_volume'])
        dat.loc[(dat['match_A'] == isid) & (dat['match_B'] == isidB),
                'inj_corr'] = spearmanr

dat.to_csv(os.path.join(path, 'match_correlations_by_source_NPV_ipsi_with_inj_corr.csv'), 
           index  = False)