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

import platform
if platform.system() == 'Windows':
    basepath = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN'
elif platform.system() == 'Darwin':
    basepath = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN'

datapath = os.path.join(basepath, 'data_files')
savepath = os.path.join(basepath, '_new_figures', 'Figure_4')

dat = pd.read_csv(os.path.join(savepath, 'TD_data_Fig4.csv'))

mod_dict = {
    'Prefrontal': ['ACAd', 'ACAv', 'ORBl','ORBm', 'ORBvl', 'ILA', 'PL', 'FRP'],
    'Medial': ['VISam', 'RSPagl', 'RSPv', 'VISpm','VISa', 'RSPd'],
    'Somatomotor': ['SSp-tr', 'SSp-ll', 'SSp-bfd', 'SSp-un',
                                                    'SSp-ul', 'MOs', 'MOp', 'SSs', 'SSp-n',
                                                    'SSp-m'],
    'Visual': ['VISal', 'VISl', 'VISli', 'VISp', 'VISpl',
                                                    'VISpor', 'VISrl'],
    'Auditory': ['AUDpo', 'AUDd', 'AUDp', 'AUDv'],
    'Lateral': ['TEa', 'PERI', 'ECT', 'GU', 'AId', 'AIv',
                                                    'AIp', 'VISC']}

source_modules = []
proj_fracs = {}
for module in mod_dict.keys():
    proj_fracs[module] = []
for isid in dat['id'].unique():
    print(isid)
    source_mod = dat[dat['id'] == isid]['module'].values[0]
    source_modules.append(source_mod)
    unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                               structure_ids = iso_desc,
                                               hemisphere_ids = [3],
                                               is_injection = False)

    for module in mod_dict.keys():
        in_targets = [ia_map[target] for target in mod_dict[module]]
        proj_frac = unionize[unionize['structure_id'].isin(in_targets)][
            'normalized_projection_volume'].sum()/unionize['normalized_projection_volume'].sum()
        proj_fracs[module].append(proj_frac)

outdat = pd.DataFrame({'image_series_id': dat['id'].unique(),
                       'source_module': source_modules})
outdat = pd.concat([outdat, pd.DataFrame(proj_fracs)], axis = 1)
outdat.to_csv(os.path.join(basepath, '_new_figures', 'Figure_4', 'all_module_projection_fractions.csv'), 
              index = False)