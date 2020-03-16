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
import platform
import json

def calculate_injection_centroid(injection_density,
                                 injection_fraction,
                                 resolution=100):
    '''
    Compute the centroid of an injection site.
    
    Parameters
    ----------
    
    injection_density: np.ndarray
        The injection density volume of an experiment

    injection_fraction: np.ndarray
        The injection fraction volume of an experiment

    '''

    # find all voxels with injection_fraction > 0
    injection_voxels = np.nonzero(injection_fraction)
    injection_density_computed = np.multiply(injection_density[injection_voxels],
                                             injection_fraction[injection_voxels]) 
    sum_density = np.sum(injection_density_computed)

    # compute centroid in CCF coordinates
    if sum_density > 0 :
        centroid = np.dot(injection_density_computed,
                          list(zip(*injection_voxels))) / sum_density * resolution
    else:
        centroid = None
    
    return centroid

def main():
    
    mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
    aapi = AnatomyApi()
    
    ss = aapi.get_summary_structure_data('id')
    structure_tree = mcc.get_structure_tree()
    ia_map = structure_tree.get_id_acronym_map()
    ai_map = {value:key for key,value in ia_map.items()}
    isocortex = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
    hipp = structure_tree.get_structures_by_acronym(['HPF'])[0]
    iso = structure_tree.descendant_ids([isocortex['id']])[0]
    hipp_desc = structure_tree.descendant_ids([hipp['id']])[0]
    ctx_ss = [strid for strid in iso if strid in ss]
    hpf_ss = [strid for strid in hipp_desc if strid in ss]
    
    td_experiments = pd.DataFrame(mcc.get_experiments(cre=['Ai75(RCL-nt)']))
    experiments =  pd.DataFrame(mcc.get_experiments())
    print(len(experiments))   
    print(args.td_id)
    if len(td_experiments[td_experiments['id'] == args.td_id]) == 0: # new data not released online
        if platform.system() == 'Windows':
            unionize_path = os.path.join(r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN_paper\alternative_unionizes',
                                         'experiment_{0}'.format(str(args.td_id)), 
                                         'output.json') #new data not online yet
        else:
            unionize_path = os.path.join(r'/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/DMN_paper/alternative_unionizes',
                                         'experiment_{0}'.format(str(args.td_id)), 
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
        inj_unionize = unionize_dat[(unionize_dat['is_injection'] == True) &
                              (unionize_dat['hemisphere_id'].isin([1,2]))]
        td_unionize = unionize_dat[(unionize_dat['is_injection'] == False) &
                                   (unionize_dat['hemisphere_id'].isin([1,2]))]
    else:
        td_unionize = mcc.get_structure_unionizes(experiment_ids = [args.td_id],
                                                is_injection = False,
                                                hemisphere_ids = [1, 2],
                                                structure_ids = ss)
        inj_unionize = mcc.get_structure_unionizes(experiment_ids = [args.td_id],
                                                is_injection = True,
                                                structure_ids = ss)
    primary = inj_unionize.sort_values(by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[0]
    
    storage_directory = aapi.get_storage_directory(args.td_id)
    data_maskA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                           'data_mask_25.nrrd'))
    injA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                     'injection_density_25.nrrd'))
    injfracA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                         'injection_fraction_25.nrrd'))
    inj_density_A = np.multiply(injA, data_maskA)
    injfracA = np.multiply(injfracA, data_maskA)
    if os.path.isfile(os.path.join(storage_directory, 'grid', 
                                   'aav_exclusion_fraction_25.nrrd')):
        exclusionA, _ = nrrd.read(os.path.join(storage_directory, 'grid', 
                                     'aav_exclusion_fraction_25.nrrd'))
        exclusion_density_A = np.multiply(exclusionA, data_maskA)
    else:
        exclusion_density_A = np.zeros_like(data_maskA)

    r_hem = inj_unionize[inj_unionize['hemisphere_id'] == 2]['projection_volume'].sum()
    l_hem = inj_unionize[inj_unionize['hemisphere_id'] == 1]['projection_volume'].sum()
    if r_hem > l_hem:
        td_experiments.loc[td_experiments['id'] == args.td_id, 'injection_z'] = 5700 - (
            td_experiments[td_experiments['id'] == args.td_id]['injection_z'] - 5700).values[0]
        inj_density_A = np.flip(inj_density_A, 2) #from R hemisphere to L
        td_unionize.loc[td_unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
        td_unionize.loc[td_unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
        td_unionize.loc[td_unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
    anchor_centroid = calculate_injection_centroid(inj_density_A, injfracA, 25)
    print(anchor_centroid)
    
    # Match data
    injB, _ = mcc.get_injection_density(experiment_id=args.wt_id)
    data_maskB, _ = mcc.get_data_mask(experiment_id=args.wt_id)
    injfracB, _ = mcc.get_injection_fraction(experiment_id=args.wt_id)
    inj_density_B = np.multiply(injB, data_maskB)
    injfracB = np.multiply(injfracB, data_maskB)
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
        injfracB = np.flip(injfracB, 2)
        match_unionize.loc[match_unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
        match_unionize.loc[match_unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
        match_unionize.loc[match_unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
    intersect = np.sum(inj_density_A[inj_density_B > 0])
    injection_coverage = np.sum(inj_density_B[inj_density_A > 0])
    if intersect > injection_coverage:
        dice = (2*float(intersect))/(float(injA.sum()) + float(injB.sum()))
    else:
        dice = (2*float(injection_coverage))/(float(injA.sum()) + float(injB.sum()))
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
    
    # threshold all at -1.5
    unionizes.loc[(np.log10(unionizes['normalized_projection_volume']) < -1.5),
                      'normalized_projection_volume'] = 0
                    
    spearmanr, _ = st.spearmanr(unionizes[unionizes['experiment_id'] == args.td_id]['normalized_projection_volume'], 
                                unionizes[unionizes['experiment_id'] == args.wt_id]['normalized_projection_volume'])
    pearsonr, _ = st.pearsonr(unionizes[unionizes['experiment_id'] == args.td_id]['normalized_projection_volume'], 
                              unionizes[unionizes['experiment_id'] == args.wt_id]['normalized_projection_volume'])
    centroid = calculate_injection_centroid(inj_density_B, injfracB, 25)
    print(centroid)
    distance = np.linalg.norm(anchor_centroid - centroid)
    print('target-defined id:', args.td_id)
    print('match id:', args.wt_id)
    print('source:', ai_map[primary])
    print('distance:', distance)
    print('injection_overlap:', overlaps)
    print('exclusion_zone_overlap:', exclusion_overlaps)
    print('fraction of match covered by td injection:', injection_fraction_covered)
    print('fraction of match covered by td exclusion zone:', exclusion_fraction_covered)
    print('dice', dice)
    print('spearman_correlation:', spearmanr)
    print('pearson_correlation:', pearsonr)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('td_id', type = int, help='target-defined experiment ID')
    parser.add_argument('wt_id', type = int, help='wild-type experiment ID')
    args = parser.parse_args()

    main()