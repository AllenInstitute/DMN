# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 13:34:58 2019

@author: jenniferwh
"""

import pandas as pd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi
import os
import json

def main():
    path = r'C:\Users\jenniferwh\Dropbox (Personal)\DMN data\correlations\core_structures'
    
    mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
    aapi = AnatomyApi()
    
    core_strs = pd.read_csv(os.path.join(path, 'ctx_hipp_proj_density_with_confidence.csv'))
    
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
    
    frac_present = []
    isids = []
    sources = []
    
    ia_map = structure_tree.get_id_acronym_map()
    ai_map = {value:key for key, value in ia_map.items()}
    for isid in experiments['id']:
        flip = False
        unionize = mcc.get_experiment_structure_unionizes(experiment_id = isid, 
                                       hemisphere_ids = [3],
                                       is_injection = True)
        inj_size = unionize[unionize['structure_id'] == 997]['projection_volume'].values[0]
        ss_dat = unionize[unionize['structure_id'].isin(ss)]
        primary = ss_dat.sort_values(by = 'projection_volume', ascending = False).reset_index()['structure_id'].values[0]
        percent_primary = ss_dat.sort_values(
                by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[0]/inj_size
        inj_unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                                is_injection = True,
                                                structure_ids = ss)
        r_hem = inj_unionize[inj_unionize['hemisphere_id'] == 2]['projection_volume'].sum()
        l_hem = inj_unionize[inj_unionize['hemisphere_id'] == 1]['projection_volume'].sum()
        if r_hem > l_hem:
            flip = True
        inj_str_unionize = mcc.get_experiment_structure_unionizes(experiment_id = isid, 
                                                                  hemisphere_ids = [3],
                                                                  is_injection = True)
        inj_size = inj_str_unionize[inj_str_unionize['structure_id'] == 997]['projection_volume'].values[0]
        ss_dat = inj_str_unionize[inj_str_unionize['structure_id'].isin(ss)]
        percent_primary = ss_dat.sort_values(
                by = 'projection_volume', ascending = False).reset_index()['projection_volume'].values[0]/inj_size
        if percent_primary > 0.5:
            unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                                    is_injection = False,
                                                    hemisphere_ids = [1, 2])
            if flip:
                unionize.loc[unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
                unionize.loc[unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
                unionize.loc[unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
            if primary in hpf_ss:
                unionize = unionize[unionize['structure_id'].isin(ctx_ss+hpf_ss)]
            else:
                unionize = unionize[unionize['structure_id'].isin(ctx_ss)]
            
            unionize.loc[(np.log10(unionize['normalized_projection_volume']) < -1.5), 'normalized_projection_volume'] = 0
            
            # Look only at core structures
            unionize['structure-hemi'] = unionize['structure_id'].map(str) + "_" + unionize['hemisphere_id'].map(str)
            core_strs['structure-hemi'] = core_strs['structure_id'].map(str) + "_" + core_strs['hemisphere_id'].map(str)
            good_strs = core_strs[(core_strs['source'] == ai_map[primary]) &
                                  (core_strs['fraction_observed'] > 0.9)]
            if len(good_strs) > 0:
                core_unionize = unionize[unionize['structure-hemi'].isin(good_strs['structure-hemi'])]
                frac_present.append(float(len(core_unionize[core_unionize['normalized_projection_volume'] > 0]))/len(core_unionize))
                isids.append(isid)
                sources.append(ai_map[primary])
        
    dat = {'image_series_id': isids, 'source': sources, 'fraction of core structures': frac_present}
    
    with open(os.path.join(path, 'core_structures_by_experiment.json'), 'w') as outfile:
        json.dump(dat, outfile, sort_keys = False, indent = 4)

if __name__ == '__main__':

    main()