# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 12:08:12 2018

@author: jenniferwh

Customizes the colors in the Allen Brain Explorer App. 

Parameters
----------
none

Returns
-------
Saves 'ontology_v2.csv' in the user-specified location. This file is used by 
the Brain Explorer App to assign colors to CCF structures.
"""

import pandas as pd
import os
from anatomy.anatomy_api import AnatomyApi
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

# Find the path to an existing ontology document that was installed through the Brain Explorer interface
path = r'C:\Users\jenniferwh\AppData\Local\Allen Institute\Brain Explorer 2\Atlases\Allen Mouse Brain Common Coordinate Framework'
dat = pd.read_csv(os.path.join(path, 'ontology_v2.csv'))
aapi = AnatomyApi()
mcc = MouseConnectivityCache()
modules = {'FRP': 'prefrontal', 
        'ACAd': 'prefrontal',
        'ACAv': 'prefrontal',
        'PL': 'prefrontal',
        'ILA': 'prefrontal',
        'ORBl': 'prefrontal',
        'ORBm': 'prefrontal',
        'ORBvl': 'prefrontal',
        'AId': 'lateral',
        'AIv': 'lateral',
        'AIp': 'lateral',
        'GU': 'lateral',
        'VISC': 'lateral',
        'TEa': 'lateral',
        'PERI': 'lateral',
        'ECT': 'lateral',
        'SSs': 'somatomotor',
        'SSp-bfd': 'somatomotor',
        'SSp-tr': 'somatomotor',
        'SSp-ll': 'somatomotor',
        'SSp-ul': 'somatomotor',
        'SSp-un': 'somatomotor',
        'SSp-n': 'somatomotor',
        'SSp-m': 'somatomotor',
        'MOp': 'somatomotor',
        'MOs': 'somatomotor',
        'VISal': 'visual',
        'VISl': 'visal',
        'VISp': 'visual',
        'VISpl': 'visual',
        'VISli': 'visual',
        'VISpor': 'visual',
        'VISrl': 'visual',
        'VISa': 'medial',
        'VISam': 'medial',
        'VISpm': 'medial',
        'RSPagl': 'medial',
        'RSPd': 'medial',
        'RSPv': 'medial',
        'AUDd': 'auditory',
        'AUDp': 'auditory',
        'AUDpo': 'auditory',
        'AUDv': 'auditory'}

def main():
    ss = aapi.get_summary_structure_data('id')
    structure_tree = mcc.get_structure_tree()
    isocortex = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
    iso = structure_tree.descendants([isocortex['id']])[0]
    ctx_ss = [struct['acronym'] for struct in iso if struct in ss]
    for structure in ctx_ss:
        if modules[structure] == 'prefrontal':
            dat.loc[dat['abbreviation'] == structure, 'red'] = 255
            dat.loc[dat['abbreviation'] == structure, 'green'] = 0
            dat.loc[dat['abbreviation'] == structure, 'blue'] = 0
        if modules[structure] == 'lateral':
            dat.loc[dat['abbreviation'] == structure, 'red'] = 255
            dat.loc[dat['abbreviation'] == structure, 'green'] = 255
            dat.loc[dat['abbreviation'] == structure, 'blue'] = 102
        if modules[structure] == 'somatomotor':
            dat.loc[dat['abbreviation'] == structure, 'red'] = 249
            dat.loc[dat['abbreviation'] == structure, 'green'] = 146
            dat.loc[dat['abbreviation'] == structure, 'blue'] = 43
        if modules[structure] == 'visual':
            dat.loc[dat['abbreviation'] == structure, 'red'] = 144
            dat.loc[dat['abbreviation'] == structure, 'green'] = 191
            dat.loc[dat['abbreviation'] == structure, 'blue'] = 249
        if modules[structure] == 'medial':
            dat.loc[dat['abbreviation'] == structure, 'red'] = 82
            dat.loc[dat['abbreviation'] == structure, 'green'] = 82
            dat.loc[dat['abbreviation'] == structure, 'blue'] = 169
        if modules[structure] == 'auditory':
            dat.loc[dat['abbreviation'] == structure, 'red'] = 124
            dat.loc[dat['abbreviation'] == structure, 'green'] = 66
            dat.loc[dat['abbreviation'] == structure, 'blue'] = 155
    dat.loc[dat['abbreviation'] == 'root', 'parent'] = 0.
    dat['parent'] = [int(value) for value in dat['parent']]
    
    dat.to_csv(os.path.join(path, 'ontology_v2.csv'), index=False)
    
if __name__ == '__main__':

    main()