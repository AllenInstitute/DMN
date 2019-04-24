#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 12:48:09 2019

@author: jenniferwh
"""
import pandas as pd
import os
import numpy as np
import scipy.stats

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi

import platform

if platform.system() == 'Darwin':
    path = '/Users/jenniferwh/Dropbox/DMN data/correlations'
elif platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox\DMN data\correlations'
    
alldat = pd.read_csv(os.path.join(path, 'matched_td_wt_cre_all.csv'))
c_by_source = pd.read_csv(os.path.join(path, 'match_correlations_by_source.csv'))

# remove duplicates
c_by_source['index_original'] = c_by_source.groupby(['match_A', 'match_B']).match_A.transform('idxmin')    
c_by_source = c_by_source[~c_by_source.duplicated(subset=['match_A', 'match_B'], keep='first')]
for isid in c_by_source['match_A'].unique():
    Bmatches = c_by_source[c_by_source['match_A'] == isid]['match_B'].values
    Amatches = c_by_source[c_by_source['match_B'] == isid]['match_A'].values
    duplicates = [match for match in Amatches if match in Bmatches]
    if len(duplicates) > 0:
        print(Bmatches)
        print(Amatches)
fail_expts = [114008926, 120280939, 180073473, 180403712, 180601025, 183174303, 183329222,
              249396394, 296047806, 299446445, 301060890, 303784745, 480069939, 482578964, 
              506947040, 514333422, 525796603, 545428296, 559878074, 638314843, 182888003] 
alldat = alldat[~alldat['match_id'].isin(fail_expts)]
c_by_source = c_by_source[~c_by_source['match_A'].isin(fail_expts)]
c_by_source = c_by_source[~c_by_source['match_B'].isin(fail_expts)]

print(len(alldat))
print(len(c_by_source))

aapi = AnatomyApi()
data = pd.read_csv(os.path.join(path, 'cla_correlations.csv'))

mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
ss_ids = aapi.get_summary_structure_data('id')

all_isids = np.unique(np.concatenate((c_by_source['match_A'].unique(), 
                                      c_by_source['match_B'].unique())))
inj_volumes = []
sum_projections = []
for experiment in all_isids: 
    unionize = mcc.get_experiment_structure_unionizes(experiment_id = experiment, 
                                           hemisphere_ids = [3])
    inj_volumes.append(unionize[(unionize['structure_id'] == 997) & 
                         (unionize['is_injection'] == True)]['projection_volume'].values[0])
    proj_unionize = unionize[(unionize['is_injection'] == False) & 
                             (unionize['structure_id'].isin(ss_ids))]
    proj_unionize.loc[(np.log10(proj_unionize['normalized_projection_volume']) < -1.5),
                       'normalized_projection_volume'] = 0
    sum_projections.append(proj_unionize['normalized_projection_volume'].sum())

X = inj_volumes
y = sum_projections
m, b, r, p, e = scipy.stats.linregress(X, y)

from scipy.optimize import curve_fit
def func_log(x, c, c2):
    return c + np.log10(x) * c2

target_func = func_log
X = inj_volumes
y = sum_projections
popt, pcov = curve_fit(target_func, X, y)

parameters = eval("scipy.stats.lognorm.fit(X)");
D, p = scipy.stats.kstest(y, "lognorm", args=parameters);
print(D, p)

for experiment in all_isids: 
    unionize = mcc.get_experiment_structure_unionizes(experiment_id = experiment, 
                                           hemisphere_ids = [3])
    inj_volumes.append(unionize[(unionize['structure_id'] == 997) & 
                         (unionize['is_injection'] == True)]['projection_volume'].values[0])
    proj_unionize = unionize[(unionize['is_injection'] == False) & 
                             (unionize['structure_id'].isin(ss_ids))]
    proj_unionize.loc[(np.log10(proj_unionize['normalized_projection_volume']) < -1.5),
                       'normalized_projection_volume'] = 0
    sum_projections.append(proj_unionize['normalized_projection_volume'].sum())