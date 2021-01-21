#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 18:07:45 2020

@author: jenniferwh
"""
import os
import pandas as pd
import platform
import statsmodels.api as sm
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi

mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                            resolution=100)
mca = MouseConnectivityApi()
structure_tree = mcc.get_structure_tree()

if platform.system() == 'Windows':
    path = r'2019 DMN'
elif platform.system() == 'Darwin':
    path = r'2019 DMN'
datpath = os.path.join(path, 'data_files')
dat = pd.read_csv(os.path.join(datpath, 
                'wt_cre_ctx_injections_DMN_and_core_projections_coefficients.csv'))


x=dat['injection dmn fraction'].values
y=dat['projection dmn fraction'].values
dmn_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(dmn_fit.summary())

x=dat['injection dmn fraction'].values
y=dat['projection core fraction'].values
core_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(core_fit.summary())

x=dat['injection dmn fraction'].values
y=dat['distance coefficient'].values
distance_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(distance_fit.summary())
