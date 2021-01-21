# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi
import scipy.stats as st
import os

path = r'DMN data/spearman_r'

mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
aapi = AnatomyApi()

ss = aapi.get_summary_structure_data('id')
structure_tree = mcc.get_structure_tree()
iso = structure_tree.descendant_ids([315])[0]
ctx_ss = [strid for strid in iso if strid in ss]

td_experiments = mcc.get_experiments(cre=['Ai75(RCL-nt)'])
td_experiments = pd.DataFrame(td_experiments)
print(len(td_experiments))

wt_expts = mcc.get_experiments(cre=False)
wt = pd.DataFrame(wt_expts)
print(len(wt))

wt['flipped_coordinates'] = [[coordinate[0], 
                             coordinate[1],
                             5700 - (coordinate[2] - 5700)] for coordinate in wt['injection-coordinates']]   

sources = []
for isid in wt['id']:
    ijstrs = wt[wt['id'] == isid]['injection-structures'].values[0]
    sources.append([item['abbreviation'] for item in ijstrs])
wt['sources'] = sources

distances = []
nearest_wt = []
overlaps = []
correlations = []
pvalues = []
for isid in td_experiments['id']:
    # find all source injection structures
    ijstrs = td_experiments[td_experiments['id'] == isid]['injection-structures'].values[0]
    sources = [item['abbreviation'] for item in ijstrs]
    # find wild type experiments with the same source injection
    match = [set(source).intersection(sources) for source in wt['sources']]
    boolmatch = [len(subset) > 0 for subset in match]
    matched_wt = wt[boolmatch]
    nearest_wt_id = 0
    max_overlap = 0
    anchor_centroid = np.array(
            td_experiments[td_experiments['id'] == isid]['injection-coordinates'].values[0])
    injA, _ = mcc.get_injection_density(experiment_id=isid)
    data_maskA, _ = mcc.get_data_mask(experiment_id=isid)
    inj_density_A = np.multiply(injA, data_maskA)
    # go through the set of wild type matches to find the one with 
    # highest overlap in injection site
    for wt_isid in matched_wt['id']:
        injB, _ = mcc.get_injection_density(experiment_id=wt_isid)
        data_maskB, _ = mcc.get_data_mask(experiment_id=wt_isid)
        inj_density_B = np.multiply(injB, data_maskB)
        inj_density_B = np.flip(inj_density_B, 2) #from R hemisphere to L
        intersect = np.sum(
            inj_density_A[inj_density_B > 0]
                          ) + np.sum(
            inj_density_B[inj_density_A > 0]
                                    )
        union = np.sum(inj_density_A + inj_density_B)
        overlap = float(intersect)/float(union)
        if overlap > max_overlap:
            max_overlap = overlap
            nearest_wt_id = wt_isid
            centroid = np.array(
                matched_wt[matched_wt['id'] == wt_isid]['flipped_coordinates'].values[0])
            dist = np.linalg.norm(anchor_centroid - centroid)
    distances.append(dist)
    nearest_wt.append(nearest_wt_id)
    unionizes = mcc.get_structure_unionizes(experiment_ids = [isid, nearest_wt_id],
                                            is_injection = False,
                                            structure_ids = ctx_ss,
                                            hemisphere_ids = [3])
    unionizes.sort_values(by='structure_id', inplace = True)
    spearmanr, pval = st.spearmanr(unionizes[unionizes['experiment_id'] == isid]['normalized_projection_volume'], 
                                   unionizes[unionizes['experiment_id'] == nearest_wt_id]['normalized_projection_volume'])
    correlations.append(spearmanr)
    pvalues.append(pval)
    
dat = pd.DataFrame({'image_series_id': isid, 'matched_wt_id': nearest_wt,
                    'distance': distances, 'overlap': overlaps, 
                    'correlation': correlations, 'pvalue': pvalues})
dat.to_csv(os.path.join(path, 'td_wild_pairs_distance_overlap_corr.csv'),
           index = False)
