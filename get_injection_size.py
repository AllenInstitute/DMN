# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:47:48 2019

@author: jenniferwh
"""

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import pandas as pd

dat = pd.read_csv(r"C:\Users\jenniferwh\Dropbox (Personal)\DMN data\correlations\all_same_primary.csv")
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
distances = {}
for isid in dat['image_series_id'].unique():
    injdat = mcc.get_structure_unionizes(experiment_ids = [isid],
                                                    is_injection = True,
                                                    structure_ids = [997],
                                                    hemisphere_ids = [3])
    dat.loc[dat['image_series_id'] == isid, 'td_injection_size']  = injdat['projection_volume'].values[0]