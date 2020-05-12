# -*- coding: utf-8 -*-
"""
Created on Wed May 22 20:09:21 2019

@author: jenniferwh
"""
import os
import numpy as np
import nrrd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

path = r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN_paper\downsampled_CAV\grid'
mcc = MouseConnectivityCache(manifest_file = '..connectivity/mouse_connectivity_manifest.json')
experiments = [521955016, 617898760, 532005897]
for n in np.arange(0, len(experiments)):
    print(n)
    for m in np.arange(n+1, len(experiments)):
        print(m)
        gridfile, _ = nrrd.read(os.path.join(path, 
                                             str(experiments[n]), 
                                             'CAV_density_25.nrrd'))
        data_mask, _ = mcc.get_data_mask(experiment_id=experiments[n])
        CAV_inj1 = np.multiply(gridfile, data_mask)
        
        expt2 = str(experiments[n+1])
        gridfile2, _ = nrrd.read(os.path.join(path, 
                                              str(experiments[m]), 
                                              'CAV_density_25.nrrd'))
        data_mask, _ = mcc.get_data_mask(experiment_id=experiments[m])
        CAV_inj2 = np.multiply(gridfile2, data_mask)
        
        # OL is always fraction of larger injection / larger injection
        if np.sum(CAV_inj1) > np.sum(CAV_inj2):
            ol = np.sum(CAV_inj1[CAV_inj2 > 0])/np.sum(CAV_inj1)
        elif np.sum(CAV_inj2) > np.sum(CAV_inj1):
            ol = np.sum(CAV_inj2[CAV_inj1 > 0])/np.sum(CAV_inj2)
        print(experiments[n])
        print(experiments[m])
        print(ol)