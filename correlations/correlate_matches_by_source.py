# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 13:34:58 2019

@author: jenniferwh
"""

import pandas as pd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi
import scipy.stats as st
import os
import json
import argparse

def main():
    path = 'cluster_code/correlations'
    
    mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
    aapi = AnatomyApi()
    
    ss = aapi.get_summary_structure_data('id')
    structure_tree = mcc.get_structure_tree()
    isocortex = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
    hipp = structure_tree.get_structures_by_acronym(['HPF'])[0]
    thal = structure_tree.get_structures_by_acronym(['TH'])[0]
    iso = structure_tree.descendant_ids([isocortex['id']])[0]
    hipp_desc = structure_tree.descendant_ids([hipp['id']])[0]
    cla = structure_tree.get_structures_by_acronym(['CLA'])[0]
    CP = structure_tree.get_structures_by_acronym(['CP'])[0]
    ctx_ss = [strid for strid in iso if strid in ss]
    hpf_ss = [strid for strid in hipp_desc if strid in ss]
    
    valid_structs = [isocortex['id'], hipp['id'], thal['id'], cla['id'], CP['id']]
    wt_expts = mcc.get_experiments(cre = False, injection_structure_ids = valid_structs)
    ctx_expts = mcc.get_experiments(cre=['Emx1-IRES-Cre', 'Rbp4-Cre_KL100'], 
                                   injection_structure_ids = [isocortex['id'],
                                                              hipp['id']])
    th_expts = mcc.get_experiments(cre=['Rbp4-Cre_KL100',
                                        'Adcyap1-2A-Cre', 'Calb1-T2A-dgCre',
                                        'Calb2-IRES-Cre', 'Cart-Tg1-Cre',
                                        'Cck-IRES-Cre', 'Chat-IRES-Cre-neo',
                                        'Crh-IRES-Cre_BL',
                                        'Efr3a-Cre_NO108', 'Foxp2-IRES-Cre',
                                        'Gal-Cre_KI97', 'Gpr26-Cre_KO250',
                                        'Grik4-Cre', 'Grm2-Cre_MR90',
                                        'Grp-Cre_KH288', 'Htr2a-Cre_KM207',
                                        'Ntrk1-IRES-Cre', 'Pdzk1ip1-Cre_KD31',
                                        'Ppp1r17-Cre_NL146', 'Prkcd-GluCla-CFP-IRES-Cre',
                                        'Scnn1a-Tg2-Cre', 'Slc17a6-IRES-Cre',
                                        'Slc17a8-IRES2-Cre', 'Slc32a1-IRES-Cre',
                                        'Tac2-IRES2-Cre', 'Vipr2-Cre_KE2',
                                        'Vipr2-IRES2-Cre'], 
                                   injection_structure_ids = [thal['id']])
    cla_expts = mcc.get_experiments(cre = ['Syt17-Cre_NO14', 'Gnb4-IRES2-Cre',
                                           'Ntng2-IRES2-Cre', 'Cux2-IRES-Cre',
                                           'Gnb4-IRES2-Cre'],
                                    injection_structure_ids = [cla['id'],
                                                               isocortex['id'],
                                                               CP['id']])
    cla_expts = [experiment for experiment in cla_expts if experiment['id'] in [
            180436360, 513773998, 296047806, 513775257, 485846989, 485903475, 187268452, 
            485902743, 513826657]]
    experiments = wt_expts + ctx_expts + th_expts + cla_expts
    experiments = pd.DataFrame(experiments)
    
    anchor_ids = []
    test_ids = []
    spearmancorr = []
    pearsoncorr = []
    distances = []
    overlaps = []
    dice = []
    absolute_overlaps = []
    absolute_overlaps_lg_sm = []
    injection_fraction_covered = []
    sum_volumes_A = []
    sum_volumes_B = []
    
    ia_map = structure_tree.get_id_acronym_map()
    potential_matches = experiments[experiments['injection_structures'].apply(lambda x: ia_map[args.source] in x)]
    potential_matches = potential_matches['id'].values
    if len(potential_matches) > 1:
        relevant_matches = experiments[experiments['id'].isin(potential_matches)]
        for isid in relevant_matches['id']:
            unionize = mcc.get_experiment_structure_unionizes(experiment_id = isid, 
                                           hemisphere_ids = [3],
                                           is_injection = True)
            inj_size = unionize[unionize['structure_id'] == 997]['projection_volume'].values[0]
            ss_dat = unionize[unionize['structure_id'].isin(ss)]
            primary = ss_dat.sort_values(by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[0]
            percent_primary = ss_dat.sort_values(
                    by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[0]/inj_size
            if primary != ia_map[args.source]:
                if len(ss_dat) > 1:
                    secondary = ss_dat.sort_values(
                            by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[1]
                    percent_secondary = ss_dat.sort_values(
                            by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[1]/inj_size
                    if secondary != ia_map[args.source]:
                        relevant_matches.drop(relevant_matches.loc[relevant_matches['id'] == isid].index, inplace = True)
                    if secondary == ia_map[args.source]:
                        if percent_primary > 0.8 or percent_secondary < 0.2:
                            relevant_matches.drop(relevant_matches.loc[relevant_matches['id'] == isid].index, inplace = True)
        matches = relevant_matches['id'].values
        relevant_matches['flipped'] = False
        for isid in relevant_matches['id']:
            inj_unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                                    is_injection = True,
                                                    structure_ids = ss)
            r_hem = inj_unionize[inj_unionize['hemisphere_id'] == 2]['projection_volume'].sum()
            l_hem = inj_unionize[inj_unionize['hemisphere_id'] == 1]['projection_volume'].sum()
            if r_hem > l_hem:
                relevant_matches.loc[relevant_matches['id'] == isid, 'flipped'] = True
                relevant_matches.loc[relevant_matches['id'] == isid, 'injection_z'] = 5700 - (
                        relevant_matches[relevant_matches['id'] == isid]['injection_z'] - 5700).values[0]
        for n in range(len(matches)):
            inj_str_unionize = mcc.get_experiment_structure_unionizes(experiment_id = matches[n], 
                                                                      hemisphere_ids = [3],
                                                                      is_injection = True)
            inj_size = inj_str_unionize[inj_str_unionize['structure_id'] == 997]['projection_volume'].values[0]
            ss_dat = inj_str_unionize[inj_str_unionize['structure_id'].isin(ss)]
            percent_primary = ss_dat.sort_values(
                    by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[0]/inj_size
            if percent_primary > 0.5:
                data_maskA, _ = mcc.get_data_mask(experiment_id=matches[n])
                injA, _ = mcc.get_injection_density(experiment_id=matches[n])
                inj_density_A = np.multiply(injA, data_maskA)
                anchor_centroid = np.array([
                    relevant_matches[relevant_matches['id'] == matches[n]]['injection_x'].values[0],
                    relevant_matches[relevant_matches['id'] == matches[n]]['injection_y'].values[0],
                    relevant_matches[relevant_matches['id'] == matches[n]]['injection_z'].values[0]])
                anchor_unionize = mcc.get_structure_unionizes(experiment_ids = [matches[n]],
                                                        is_injection = False,
                                                        hemisphere_ids = [1, 2])
                if ia_map[args.source] in hpf_ss:
                    anchor_unionize = anchor_unionize[anchor_unionize['structure_id'].isin(ctx_ss+hpf_ss)]
                else:
                    anchor_unionize = anchor_unionize[anchor_unionize['structure_id'].isin(ctx_ss)]
                if relevant_matches[relevant_matches['id'] == matches[n]]['flipped'].values[0] == True:
                    inj_density_A = np.flip(inj_density_A, 2) #from R hemisphere to L
                    anchor_unionize.loc[anchor_unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
                    anchor_unionize.loc[anchor_unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
                    anchor_unionize.loc[anchor_unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
                
                for m in range(len(matches) - n - 1):
                    ix = m+n+1
                    anchor_ids.append(int(matches[n]))
                    test_ids.append(int(matches[ix]))
                    injB, _ = mcc.get_injection_density(experiment_id=matches[ix])
                    data_maskB, _ = mcc.get_data_mask(experiment_id=matches[ix])
                    inj_density_B = np.multiply(injB, data_maskB)
                    test_inj_str_unionize = mcc.get_experiment_structure_unionizes(experiment_id = matches[ix], 
                                                                                   hemisphere_ids = [3],
                                                                                   is_injection = True)
                    test_inj_size = test_inj_str_unionize[test_inj_str_unionize['structure_id'] == 997]['projection_volume'].values[0]
                    test_unionize = mcc.get_structure_unionizes(experiment_ids = [matches[ix]],
                                                            is_injection = False,
                                                            hemisphere_ids = [1, 2])
                    if ia_map[args.source] in hpf_ss:
                        test_unionize = test_unionize[test_unionize['structure_id'].isin(ctx_ss+hpf_ss)]
                    else:
                        test_unionize = test_unionize[test_unionize['structure_id'].isin(ctx_ss)]
                    if relevant_matches[relevant_matches['id'] == matches[ix]]['flipped'].values[0] == True:
                        inj_density_B = np.flip(inj_density_B, 2) #from R hemisphere to L
                        test_unionize.loc[test_unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
                        test_unionize.loc[test_unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
                        test_unionize.loc[test_unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
                    
                    intersect = np.sum(inj_density_A[inj_density_B > 0])
                    injection_coverage = np.sum(inj_density_B[inj_density_A > 0])
                    
                    if intersect > 0:
                        absolute_overlaps.append(intersect)
                        absolute_overlaps_lg_sm.append(injection_coverage)
                        overlaps.append(float(intersect)/float(inj_density_A.sum()))
                        injection_fraction_covered.append(float(injection_coverage)/float(inj_density_B.sum()))
                    else:
                        absolute_overlaps.append(0)
                        absolute_overlaps_lg_sm.append(0)
                        overlaps.append(0)
                        injection_fraction_covered.append(0)
                    if intersect > injection_coverage:
                        dice.append((2*float(intersect))/(float(injA.sum()) + float(injB.sum())))
                    else:
                        dice.append((2*float(injection_coverage))/(float(injA.sum()) + float(injB.sum())))    
                    unionizes = pd.concat([anchor_unionize, test_unionize])
                    unionizes.sort_values(by=['hemisphere_id', 'structure_id'], inplace = True)
                    '''
                    # threshold wt at log10(weight) = -3.5, cre at log10(weight) = -1.5
                    if relevant_matches[relevant_matches['id'] == matches[n]]['transgenic_line'].values[0] is None:
                        unionizes.loc[(unionizes['experiment_id'] == matches[n]) &
                                  (np.log10(unionizes['normalized_projection_volume']) < -3.5),
                                  'normalized_projection_volume'] = 0
                    else:
                        unionizes.loc[(unionizes['experiment_id'] == matches[n]) &
                                  (np.log10(unionizes['normalized_projection_volume']) < -1.5),
                                  'normalized_projection_volume'] = 0
                    '''
                    # threshold all unionizes at -1.5
                    unionizes.loc[(np.log10(unionizes['normalized_projection_volume']) < -1.5),
                                  'normalized_projection_volume'] = 0
                                       
                    test_centroid = np.array([
                        relevant_matches[relevant_matches['id'] == matches[ix]]['injection_x'].values[0],
                        relevant_matches[relevant_matches['id'] == matches[ix]]['injection_y'].values[0],
                        relevant_matches[relevant_matches['id'] == matches[ix]]['injection_z'].values[0]])
                         
                    distances.append(np.linalg.norm(anchor_centroid - test_centroid))
                    spearmanr, _ = st.spearmanr(unionizes[unionizes['experiment_id'] == matches[n]]['normalized_projection_volume'], 
                                                unionizes[unionizes['experiment_id'] == matches[ix]]['normalized_projection_volume'])
                    pearsonr, _ = st.pearsonr(unionizes[unionizes['experiment_id'] == matches[n]]['normalized_projection_volume'], 
                                              unionizes[unionizes['experiment_id'] == matches[ix]]['normalized_projection_volume'])
                    spearmancorr.append(spearmanr)
                    pearsoncorr.append(pearsonr)
                    sum_volumes_A.append(unionizes[unionizes['experiment_id'] == matches[n]]['normalized_projection_volume'].sum())
                    sum_volumes_B.append(unionizes[unionizes['experiment_id'] == matches[ix]]['normalized_projection_volume'].sum())
            
    dat = {'match_A': anchor_ids, 'match_B': test_ids,
           'distance': distances, 'absolute_overlap': absolute_overlaps,
           'absolute_overlap_lg_over_small': absolute_overlaps_lg_sm,
           'injection_overlap': overlaps, 
           'fraction of match covered by td injection': injection_fraction_covered,
           'dice_coefficient': dice,
           'spearman_correlation': spearmancorr, 
           'pearson_correlation': pearsoncorr, 'sum_NPV_A': sum_volumes_A,
           'sum_NPV_B': sum_volumes_B, 'source': args.source}
    
    with open(os.path.join(path, 'output', 'matches_by_source', 'NPV', '{}.json'.format(args.source)), 'w') as outfile:
        json.dump(dat, outfile, sort_keys = False, indent = 4)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('source', type = str, help='source')
    args = parser.parse_args()

    main()
