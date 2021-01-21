#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 20:03:40 2020

@author: jenniferwh
"""

import pandas as pd
import numpy as np
import os
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi
import json

aapi = AnatomyApi()
ss = aapi.get_summary_structure_data('id')
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
structure_tree = mcc.get_structure_tree()
isocortex = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
iso = structure_tree.descendant_ids([isocortex['id']])[0]
iso = [structure for structure in iso if structure in ss]
ia_map = structure_tree.get_id_acronym_map()
ai_map = {value:key for key, value in ia_map.items()}
ctx_strs = [ai_map[structure] for structure in iso]

import platform
if platform.system() == 'Darwin':
    basepath = '2019 DMN'
elif platform.system() == 'Windows':
    basepath = ''
datpath = os.path.join(basepath, 'data_files')
savepath = os.path.join(basepath, '_new_figures', 'Figure_4')

td_dataset = pd.read_csv(os.path.join(basepath, 'target_defined_dataset.csv'))
td_dataset = td_dataset[td_dataset['include'] == 'yes']

print(len(td_dataset))
print(len(td_dataset['source'].unique()))

c_by_source = pd.read_csv(os.path.join(datpath, 'good_wt_correlations.csv'))
print(len(c_by_source))
alldat = pd.read_csv(os.path.join(datpath, 'good_td_wt_correlations.csv'))
print(len(alldat))
td_dat = pd.read_csv(os.path.join(datpath, 'good_td_td_correlations.csv'))
print(len(td_dat))

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)

groups = []
for source in td_dat['source'].unique():
    matches = []
    wt_sets = []
    counter = 0
    valid_td = []
    if len(td_dat[(td_dat['source'] == source) &
                    (td_dat['experiment_type'] == 'in-out')]) > 0:
        print(source)
        # find sets of in-out targets in this source
        td_matches = td_dat[(td_dat['source'] == source) & 
                            (td_dat['experiment_type'] == 'in-out')]
        # find wt matches for this in-out pair
        for ix, row in td_matches.iterrows():
            wtA = alldat[alldat['image_series_id'] == row['image_series_id']][
                'match_id'].unique() #wt matches to TD1
            wtB = alldat[alldat['image_series_id'] == row['match_id']][
                'match_id'].unique() #wt matches to TD2
            wt_matches = [isid for isid in wtA if isid in wtB] #wt matches for both
            if len(wt_matches) < 1:
                print('no wt matches found')
            else:
                counter += 1
            # identify in and out experiments
            if row['td_CAV_percent_DMN'] > 50:
                in_exp = row['image_series_id']
                out_exp = row['match_id']
                if row['match_CAV_percent_DMN'] > 50:
                    print('something is wrong with in-out by target detection')
            else:
                in_exp = row['match_id']
                out_exp = row['image_series_id']
                if row['match_CAV_percent_DMN'] < 50:
                    print('something is wrong with in-out by target detection')
            valid_td = np.append(valid_td, [in_exp])
            valid_td = np.append(valid_td, [out_exp])
            
            # find all TD matches for the "in" experiment
            more_td_matches = np.unique(np.concatenate((
                td_dat[td_dat['image_series_id'] == in_exp]['match_id'].unique(),
                td_dat[td_dat['match_id'] == in_exp]['image_series_id'].unique())))
            # only keep additional TD matches if they match the out expt and at least one WT
            for match in more_td_matches:
                wt = alldat[alldat['image_series_id'] == match][
                    'match_id'].unique() 
                if len([match for match in wt if match in wt_matches]) > 0:
                    tdm = td_dat[td_dat['image_series_id'] == match][
                        'match_id'].unique()
                    if out_exp not in tdm:
                        tdm = td_dat[td_dat['match_id'] == match][
                            'image_series_id'].unique()
                        if out_exp in tdm:
                            valid_td = np.append(valid_td, [match])
                    elif out_exp in tdm:
                        valid_td = np.append(valid_td, [match])
                 
            # find all TD matches for the "out" experiment
            more_td_matches = np.unique(np.concatenate((
                td_dat[td_dat['image_series_id'] == out_exp]['match_id'].unique(),
                td_dat[td_dat['match_id'] == out_exp]['image_series_id'].unique())))
            # only keep TD matches if they match the out expt and at least one WT
            for match in more_td_matches:
                wt = alldat[alldat['image_series_id'] == match][
                    'match_id'].unique()
                if len([match for match in wt if match in wt_matches]) > 0:
                    tdm = td_dat[td_dat['image_series_id'] == match][
                        'match_id'].unique()
                    if out_exp not in tdm:
                        tdm = td_dat[td_dat['match_id'] == match][
                            'image_series_id'].unique()
                        if out_exp in tdm:
                            valid_td = np.append(valid_td, [in_exp])
                    elif out_exp in tdm:
                        valid_td = np.append(valid_td, [in_exp])
            
            valid_td = np.unique(valid_td)
            for td in valid_td:
                wt = alldat[alldat['image_series_id'] == match][
                    'match_id'].unique()
                wt_matches = [int(isid) for isid in wt_matches if isid in wt] #remove wt matches
                # that do not match with additional td experiment
            wt_sets.append(wt_matches)
            matches.append([int(value) for value in valid_td])
            counter += 1
        
        matches = [list(match) for match in matches]
        wt_sets = [list(wts) for wts in wt_sets]
        groups.append({'source': source, 'td_sets': np.unique(matches), 
                       'wt_sets': np.unique(wt_sets)})
        print('Num groups', len(np.unique(matches)))
        print('Num wt experiments in group', len(np.unique(wt_matches)))
with open(os.path.join(savepath, 'matched_sets.json'), 'w') as outfile:
    json.dump(groups, outfile, sort_keys = False, indent = 4, cls = MyEncoder)
    outfile.write("\n")  # Add newline cause Py JSON does not
        
