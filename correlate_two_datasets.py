# -*- coding: utf-8 -*-
"""
Created on Thu Nov 01 14:54:09 2018

@author: jenniferwh
"""

import pandas as pd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi
import scipy.stats as st
import os
import argparse
import nrrd

def main():
    
    mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
    aapi = AnatomyApi()
    
    ss = aapi.get_summary_structure_data('id')
    structure_tree = mcc.get_structure_tree()
    isocortex = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
    hipp = structure_tree.get_structures_by_acronym(['HPF'])[0]
    iso = structure_tree.descendant_ids([isocortex['id']])[0]
    hipp_desc = structure_tree.descendant_ids([hipp['id']])[0]
    ctx_ss = [strid for strid in iso if strid in ss]
    hpf_ss = [strid for strid in hipp_desc if strid in ss]
    
    td_experiments = pd.DataFrame(mcc.get_experiments(cre=['Ai75(RCL-nt)']))
    print(len(td_experiments))
    experiments =  pd.DataFrame(mcc.get_experiments())
    print(len(experiments))   
    unionize = mcc.get_experiment_structure_unionizes(experiment_id = args.td_id, 
                                           hemisphere_ids = [3],
                                           is_injection = True)
    ss_dat = unionize[unionize['structure_id'].isin(ss)]
    primary = ss_dat.sort_values(by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[0]
    # target-defined data
    td_unionize = mcc.get_structure_unionizes(experiment_ids = [args.td_id],
                                                is_injection = False,
                                                hemisphere_ids = [1, 2])
    storage_directory = aapi.get_storage_directory(args.td_id)
    data_maskA, _ = mcc.get_data_mask(experiment_id=args.td_id)
    injA, _ = mcc.get_injection_density(experiment_id=args.td_id)
    inj_density_A = np.multiply(injA, data_maskA)
    if os.path.isfile(os.path.join(storage_directory, 'grid', 
                                   'aav_exclusion_fraction_25.nrrd')):
        exclusionA, _ = nrrd.read(os.path.join(storage_directory, 'grid', 
                                     'aav_exclusion_fraction_25.nrrd'))
        exclusion_density_A = np.multiply(exclusionA, data_maskA)
    else:
        exclusion_density_A = np.zeros_like(data_maskA)
        # iterate through the set of wild type matches 
    inj_unionize = mcc.get_structure_unionizes(experiment_ids = [args.td_id],
                                                is_injection = True,
                                                structure_ids = ss)
    r_hem = inj_unionize[inj_unionize['hemisphere_id'] == 2]['projection_volume'].sum()
    l_hem = inj_unionize[inj_unionize['hemisphere_id'] == 1]['projection_volume'].sum()
    if r_hem > l_hem:
        td_experiments.loc[td_experiments['id'] == args.td_id, 'injection_z'] = 5700 - (
                    td_experiments[td_experiments['id'] == args.td_id]['injection_z'] - 5700).values[0]
        inj_density_A = np.flip(inj_density_A, 2) #from R hemisphere to L
        td_unionize.loc[td_unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
        td_unionize.loc[td_unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
        td_unionize.loc[td_unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
    anchor_centroid = np.array([
            td_experiments[td_experiments['id'] == args.td_id]['injection_x'].values[0],
            td_experiments[td_experiments['id'] == args.td_id]['injection_y'].values[0],
            td_experiments[td_experiments['id'] == args.td_id]['injection_z'].values[0]])
    
    # Match data
    injB, _ = mcc.get_injection_density(experiment_id=args.wt_id)
    data_maskB, _ = mcc.get_data_mask(experiment_id=args.wt_id)
    inj_density_B = np.multiply(injB, data_maskB)
    inj_unionize = mcc.get_structure_unionizes(experiment_ids = [args.wt_id],
                                               is_injection = True,
                                                structure_ids = ss)
    match_unionize = mcc.get_structure_unionizes(experiment_ids = [args.wt_id],
                                            is_injection = False,
                                            hemisphere_ids = [1, 2])
    r_hem = inj_unionize[inj_unionize['hemisphere_id'] == 2]['projection_volume'].sum()
    l_hem = inj_unionize[inj_unionize['hemisphere_id'] == 1]['projection_volume'].sum()
    if r_hem > l_hem:
        experiments.loc[experiments['id'] == args.wt_id, 'injection_z'] = 5700 - (
                    experiments[experiments['id'] == args.wt_id]['injection_z'] - 5700).values[0]
        inj_density_B = np.flip(inj_density_B, 2) #from R hemisphere to L
        match_unionize.loc[match_unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
        match_unionize.loc[match_unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
        match_unionize.loc[match_unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
    intersect = np.sum(inj_density_A[inj_density_B > 0])
    injection_coverage = np.sum(inj_density_B[inj_density_A > 0])
    if exclusion_density_A.sum() > 0:
        exclusion_intersect = np.sum(exclusion_density_A[inj_density_B > 0])
        exclusion_coverage = np.sum(inj_density_B[exclusion_density_A > 0])
    else:
        exclusion_intersect = 0
        exclusion_coverage = 0
        exclusion_density_A = np.zeros_like(data_maskA)
    if intersect > 0:
        overlaps = float(intersect)/float(inj_density_A.sum())
        injection_fraction_covered = float(injection_coverage)/float(inj_density_B.sum())
    else:
        overlaps = 0
        injection_fraction_covered = 0
    if exclusion_intersect > 0:
        exclusion_overlaps = float(exclusion_intersect)/float(exclusion_density_A.sum())
        exclusion_fraction_covered = float(exclusion_coverage)/float(inj_density_B.sum())
    else:
        exclusion_overlaps = 0
        exclusion_fraction_covered = 0

    unionizes = pd.concat([td_unionize, match_unionize])
    unionizes.sort_values(by=['hemisphere_id', 'structure_id'], inplace = True)
    if primary in hpf_ss:
        unionizes = unionizes[unionizes['structure_id'].isin(ctx_ss+hpf_ss)]
    else:
        unionizes = unionizes[unionizes['structure_id'].isin(ctx_ss)]
    
    # threshold wt and cre at -1.5, cre and target-defined at -2.5
    unionizes.loc[(unionizes['experiment_id'] == args.wt_id) &
                      (np.log10(unionizes['normalized_projection_volume']) < -1.5),
                      'normalized_projection_volume'] = 0
    unionizes.loc[(unionizes['experiment_id'] == args.td_id) &
                      (np.log10(unionizes['normalized_projection_volume']) < -2.5),
                      'normalized_projection_volume'] = 0
                    
    spearmanr, _ = st.spearmanr(unionizes[unionizes['experiment_id'] == args.td_id]['normalized_projection_volume'], 
                                unionizes[unionizes['experiment_id'] == args.wt_id]['normalized_projection_volume'])
    pearsonr, _ = st.pearsonr(unionizes[unionizes['experiment_id'] == args.td_id]['normalized_projection_volume'], 
                              unionizes[unionizes['experiment_id'] == args.wt_id]['normalized_projection_volume'])
    centroid = np.array([
            experiments[experiments['id'] == args.wt_id]['injection_x'].values[0],
            experiments[experiments['id'] == args.wt_id]['injection_y'].values[0],
            experiments[experiments['id'] == args.wt_id]['injection_z'].values[0]])
    distance = np.linalg.norm(anchor_centroid - centroid)
    print('target-defined id:', args.td_id)
    print('match id:', args.wt_id)
    print('source:', td_experiments[td_experiments['id'] == args.td_id]['structure_abbrev'].values[0])
    print('distance:', distance)
    print('injection_overlap:', overlaps)
    print('exclusion_zone_overlap:', exclusion_overlaps)
    print('fraction of match covered by td injection:', injection_fraction_covered)
    print('fraction of match covered by td exclusion zone:', exclusion_fraction_covered)
    print('spearman_correlation:', spearmanr)
    print('pearson_correlation:', pearsonr)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('td_id', type = int, help='target-defined experiment ID')
    parser.add_argument('wt_id', type = int, help='wild-type experiment ID')
    args = parser.parse_args()

    main()