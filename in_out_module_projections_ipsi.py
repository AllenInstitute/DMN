# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 19:37:48 2020

@author: jenniferwh
"""

import os
import pandas as pd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi
aapi = AnatomyApi()
ss = aapi.get_summary_structure_data('id')
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
iso = st.get_structures_by_acronym(['Isocortex'])[0]
iso_desc = st.descendants([iso['id']])[0]
iso_desc = [structure['id'] for structure in iso_desc]
iso_desc = [structure for structure in iso_desc if structure in ss]

def check_for_flip(unionize):
    r_hem = unionize[unionize['hemisphere_id'] == 2]['projection_volume'].sum()
    l_hem = unionize[unionize['hemisphere_id'] == 1]['projection_volume'].sum()
    if r_hem > l_hem:
        return 2
    else:
        return 1

import platform
if platform.system() == 'Windows':
    basepath = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN'
elif platform.system() == 'Darwin':
    basepath = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN'

datapath = os.path.join(basepath, 'data_files')
savepath = os.path.join(basepath, '_new_figures', 'Figure_3')

dat = pd.read_csv(os.path.join(savepath, 'curated_cre_layer_data.csv'))

mod_dict = {
    'PFC': ['ACAd', 'ACAv', 'ORBl','ORBm', 'ORBvl', 'ILA', 'PL', 'FRP'],
    'Medial': ['VISam', 'RSPagl', 'RSPv', 'VISpm','VISa', 'RSPd'],
    'Somatomotor': ['SSp-tr', 'SSp-ll', 'SSp-bfd', 'SSp-un',
                                                    'SSp-ul', 'MOs', 'MOp', 'SSs', 'SSp-n',
                                                    'SSp-m'],
    'Visual': ['VISal', 'VISl', 'VISli', 'VISp', 'VISpl',
                                                    'VISpor', 'VISrl'],
    'Auditory': ['AUDpo', 'AUDd', 'AUDp', 'AUDv'],
    'Lateral': ['TEa', 'PERI', 'ECT', 'GU', 'AId', 'AIv',
                                                    'AIp', 'VISC']}

sources = []
modules = []
proj_fracs = []
for isid in dat['image_series_id'].unique():
    print(isid)
    source = dat[dat['image_series_id'] == isid]['Exp Source'].values[0]
    print(source)
    sources.append(source)
    module = dat[dat['image_series_id'] == isid]['Exp Module'].values[0]
    print(module)
    modules.append(module)
    in_targets = [ia_map[target] for target in mod_dict[module]]
    inj_unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                           structure_ids = iso_desc,
                                           hemisphere_ids = [1,2],
                                           is_injection = True)
    unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                           structure_ids = iso_desc,
                                           hemisphere_ids = [check_for_flip(inj_unionize)],
                                           is_injection = False)
    proj_frac = unionize[unionize['structure_id'].isin(in_targets)][
        'normalized_projection_volume'].sum()/unionize['normalized_projection_volume'].sum()
    proj_fracs.append(proj_frac)

outdat = pd.DataFrame({'image_series_id': dat['image_series_id'].unique(),
                       'source': sources,
                       'module_proj_frac': proj_fracs})
outdat.to_csv(os.path.join(basepath, '_new_figures', 'Figure_3', 'module_projection_fractions_ipsi.csv'), index = False)