#!/shared/utils.x86_64/python-2.7/bin/python
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
import json
import argparse
import nrrd
import platform

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
    isocortex = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
    hipp = structure_tree.get_structures_by_acronym(['HPF'])[0]
    thal = structure_tree.get_structures_by_acronym(['TH'])[0]
    cla = structure_tree.get_structures_by_acronym(['CLA'])[0]
    CP = structure_tree.get_structures_by_acronym(['CP'])[0]
    iso = structure_tree.descendant_ids([isocortex['id']])[0]
    hipp_desc = structure_tree.descendant_ids([hipp['id']])[0]
    ctx_ss = [strid for strid in iso if strid in ss]
    hpf_ss = [strid for strid in hipp_desc if strid in ss]

    ia_map = structure_tree.get_id_acronym_map()
    ai_map = {value:key for key, value in ia_map.items()}
    
    td_experiments = mcc.get_experiments(cre=['Ai75(RCL-nt)'])
    td_experiments = pd.DataFrame(td_experiments)
    print(len(td_experiments))
    wt_expts = mcc.get_experiments(cre = False, injection_structure_ids = [isocortex['id'],
                                                                           hipp['id'],
                                                                           cla['id']])
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

    # find target-defined experiment source
    if len(td_experiments[td_experiments['id'] == args.id]) == 0: # new data not released online
        if platform.system() == 'Windows':
             unionize_path = os.path.join(r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN_paper\alternative_unionizes',
                                         'experiment_{0}'.format(str(args.id)), 
                                         'output.json') #new data not online yet
        else:
            unionize_path = os.path.join(r'/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/DMN_paper/alternative_unionizes',
                                         'experiment_{0}'.format(str(args.id)), 
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
        injdat = unionize_dat[(unionize_dat['is_injection'] == True) &
                              (unionize_dat['structure_id'].isin(ss)) &
                              (unionize_dat['hemisphere_id'] == 3)]
        sources = injdat['structure_id'].unique()
        print([ai_map[source] for source in sources])
        source_abbrev = ai_map[injdat[injdat['projection_volume'] == injdat['projection_volume'].max()]['structure_id'].values[0]]
        if source_abbrev == 'fiber tracts':
            source_abbrev = ai_map[injdat[injdat['projection_volume'] == injdat['projection_volume'].max()]['structure_id'].values[1]]
        unionize = unionize_dat[(unionize_dat['is_injection'] == True) &
                              (unionize_dat['hemisphere_id'] == 3)]
        td_unionize = unionize_dat[(unionize_dat['is_injection'] == False) &
                                   (unionize_dat['hemisphere_id'] == 1)]
    else: #most experiments in sdk
        sources = td_experiments[td_experiments['id'] == args.id]['injection_structures'].values[0]
        source_abbrev = td_experiments[td_experiments['id'] == args.id]['structure_abbrev'].values[0]
        if source_abbrev == 'fiber tracts':
            source_abbrev = ai_map[injdat[injdat['projection_volume'] == injdat['projection_volume'].max()]['structure_id'].values[1]]
        unionize = mcc.get_experiment_structure_unionizes(experiment_id = args.id,
                                       hemisphere_ids = [3],
                                       is_injection = True)
        td_unionize = mcc.get_structure_unionizes(experiment_ids = [args.id],
                                            is_injection = False,
                                            hemisphere_ids = [1])
    td_unionize = td_unionize[['experiment_id', 'hemisphere_id', 'is_injection', 
                               'structure_id', 'projection_volume', 'normalized_projection_volume']]
    # For some reason a few sources are missing from alternative-unioniize file. Hack to make it work
    print(len(td_unionize))
    '''
    td_unionize = td_unionize.append({'experiment_id': args.id,
                        'hemisphere_id': 2,
                        'is_injection': False,
                        'structure_id': 44,
                        'projection_volume': 0,
                        'normalized_projection_volume': 0},
        ignore_index = True)

    td_unionize = td_unionize.append({'experiment_id': args.id,
                        'hemisphere_id': 2,
                        'is_injection': False,
                        'structure_id': 312782574,
                        'projection_volume': 0,
                        'normalized_projection_volume': 0},
        ignore_index = True)
    print(len(td_unionize))'''
    cla_sources = [104, 119, 583, 672]
    if len(set(sources).intersection(cla_sources)) > 0: #claustrum injection
        matches = pd.DataFrame(cla_expts)
    else:
        match = [set(source_ids).intersection(sources) for source_ids in experiments['injection_structures']]
        boolmatch = [len(subset) > 0 for subset in match]
        matches = experiments[boolmatch]
        print(matches['id'].unique())
    
    ss_dat = unionize[unionize['structure_id'].isin(ss)]
    primary = ss_dat.sort_values(by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[0]
    
    # get grid data for target-defined experiment
    storage_directory = aapi.get_storage_directory(args.id)
    data_maskA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                           'data_mask_25.nrrd'))
    injA, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                     'injection_density_25.nrrd'))
    INJFRAC, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                        'injection_fraction_25.nrrd'))
    anchor_centroid = calculate_injection_centroid(injA, INJFRAC, 25)
    inj_density_A = np.multiply(injA, data_maskA)
    if os.path.isfile(os.path.join(storage_directory, 'grid', 
                                   'aav_exclusion_fraction_25.nrrd')):
        exclusionA, _ = nrrd.read(os.path.join(storage_directory, 'grid', 
                                     'aav_exclusion_fraction_25.nrrd'))
        exclusion_density_A = np.multiply(exclusionA, data_maskA)
    else:
        exclusion_density_A = np.zeros_like(data_maskA)
        
    # flip right hemisphere injections to the left side
    matches['flipped'] = False
    for isid in matches['id']:
        inj_unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                                    is_injection = True,
                                                    structure_ids = ss)
        r_hem = inj_unionize[inj_unionize['hemisphere_id'] == 2]['projection_volume'].sum()
        l_hem = inj_unionize[inj_unionize['hemisphere_id'] == 1]['projection_volume'].sum()
        if r_hem > l_hem:
            matches.loc[matches['id'] == isid, 'flipped'] = True

    spearmancorr = []
    pearsoncorr = []
    distances = []
    overlaps = []
    exclusion_overlaps = []
    injection_fraction_covered = []
    exclusion_fraction_covered = []
    match_isids = []
    
    # iterate through the set of wild type matches and calculate correlations
    for isid in matches['id']:
        injB, _ = mcc.get_injection_density(experiment_id=isid)
        inj_frac_B, _ = mcc.get_injection_fraction(experiment_id=isid)
        data_maskB, _ = mcc.get_data_mask(experiment_id=isid)
        inj_density_B = np.multiply(injB, data_maskB)
        inj_unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                                   hemisphere_ids = [3],
                                                    is_injection = True)
        match_unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                                is_injection = False,
                                                hemisphere_ids = [1, 2])
        # flip unionize data for right side injections
        if matches[matches['id'] == isid]['flipped'].values[0] == True:
            injB = np.flip(injB, 2) #from R hemisphere to L
            inj_density_B = np.flip(inj_density_B, 2) #from R hemisphere to L
            inj_frac_B = np.flip(inj_frac_B, 2) #from R hemisphere to L
            match_unionize.loc[match_unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
            match_unionize.loc[match_unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
            match_unionize.loc[match_unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
        match_unionize = match_unionize[['experiment_id', 'hemisphere_id', 'is_injection',
                               'structure_id', 'projection_volume', 'normalized_projection_volume']]
        match_unionize = match_unionize[match_unionize['hemisphere_id'] == 1]
        intersect = np.sum(inj_density_A[inj_density_B > 0])
        injection_coverage = np.sum(inj_density_B[inj_density_A > 0])
        if exclusion_density_A.sum() > 0:
            exclusion_intersect = np.sum(exclusion_density_A[inj_density_B > 0])
            exclusion_coverage = np.sum(inj_density_B[exclusion_density_A > 0])
        else:
            exclusion_intersect = 0
            exclusion_coverage = 0
            exclusion_density_A = np.zeros_like(data_maskA)
        #if intersect > 0 or exclusion_intersect > 0:
        match_isids.append(isid)
        if intersect > 0:
            overlaps.append(float(intersect)/float(inj_density_A.sum()))
            injection_fraction_covered.append(float(injection_coverage)/float(inj_density_B.sum()))
        else:
            overlaps.append(0)
            injection_fraction_covered.append(0)
        if exclusion_intersect > 0:
            exclusion_overlaps.append(float(exclusion_intersect)/float(exclusion_density_A.sum()))
            exclusion_fraction_covered.append(float(exclusion_coverage)/float(inj_density_B.sum()))
        else:
            exclusion_overlaps.append(0)
            exclusion_fraction_covered.append(0)
        unionizes = pd.concat([td_unionize, match_unionize])
        unionizes.sort_values(by=['hemisphere_id', 'structure_id'], inplace = True)
        if primary in hpf_ss:
            unionizes = unionizes[unionizes['structure_id'].isin(ctx_ss+hpf_ss)]
        else:
            unionizes = unionizes[unionizes['structure_id'].isin(ctx_ss)]

        # threshold at -1.5 for all
        unionizes.loc[np.log10(unionizes['normalized_projection_volume']) < -1.5,
                          'normalized_projection_volume'] = 0
        print(len(unionizes[unionizes['experiment_id'] == args.id]))
        print(len(unionizes[unionizes['experiment_id'] == isid]))
        print([structure for structure in unionizes[(unionizes['experiment_id'] == isid)]['structure_id'].unique() if structure not in unionizes[
                                                    (unionizes['experiment_id'] == args.id)]['structure_id'].values])
        spearmanr, _ = st.spearmanr(unionizes[unionizes['experiment_id'] == args.id]['normalized_projection_volume'],
                                    unionizes[unionizes['experiment_id'] == isid]['normalized_projection_volume'])
        pearsonr, _ = st.pearsonr(unionizes[unionizes['experiment_id'] == args.id]['normalized_projection_volume'], 
                                  unionizes[unionizes['experiment_id'] == isid]['normalized_projection_volume'])
        spearmancorr.append(spearmanr)
        pearsoncorr.append(pearsonr)
        centroid = calculate_injection_centroid(injB, inj_frac_B, 25)
        distances.append(np.linalg.norm(anchor_centroid - centroid))
        dat = {'image_series_id': args.id, 'match_id':  match_isids,
           'distance': distances, 'injection_overlap': overlaps, 
           'exclusion_zone_overlap': exclusion_overlaps,
           'fraction of match covered by td injection': injection_fraction_covered,
           'fraction of match covered by td exclusion zone': exclusion_fraction_covered,
           'spearman_correlation': spearmancorr, 
           'pearson_correlation': pearsoncorr,
           'source': source_abbrev}
    savepath = '/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/cluster_code/correlations/output/NPV_ipsi'
    if platform.system() == 'Windows':
        savepath = r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\cluster_code\correlations\output\NPV_ipsi'
    with open(os.path.join(savepath, '{}.json'.format(args.id)), 'w') as outfile:
        json.dump(dat, outfile, sort_keys = False, indent = 4)
    print(dat)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('id', type = int, help='experiment ID')
    args = parser.parse_args()

    main()