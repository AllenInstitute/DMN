# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 15:54:21 2020

@author: jenniferwh
"""

from anatomy.anatomy_api import AnatomyApi
import os
import numpy as np
import pandas as pd
import nrrd
import seaborn as sns
import matplotlib.pyplot as plt
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json',
                              resolution = 100)
st = mcc.get_structure_tree()
iso = st.get_structures_by_acronym(['Isocortex'])[0]
iso_mask = mcc.get_structure_mask(iso['id'])[0]
aapi = AnatomyApi()
basepath = r'2019 DMN'
datpath = os.path.join(basepath, 'data_files')
dat = pd.read_csv(os.path.join(datpath, 'failed_AAV_injections_for_grid_threshold.csv'))
savepath = os.path.join(basepath, '_new_figures', 'Figure_2')

grid_values = []
for isid in [476046249, 530001582]:
    gridpath = aapi.get_storage_directory(isid)
    try:
        grid, _ = nrrd.read(os.path.join(gridpath, 'grid', 'projection_density_100.nrrd'))
        data_mask, _ = nrrd.read(os.path.join(gridpath, 'grid', 'data_mask_100.nrrd'))
        grid = grid * data_mask
        grid = grid*iso_mask
        grid_values.append(grid[np.where(grid>0)].flatten())
    except:
        print('N/A')

grid_values = np.array([item for sublist in grid_values for item in sublist])
np.min(grid_values)
np.max(grid_values)
fig, ax = plt.subplots()
sns.distplot(np.log10(grid_values), kde = False, ax = ax)
ax.set_xlabel('Log$_{10}$(False Positive Value)')
ax.set_ylabel('Number of Values')
plt.show()
plt.savefig(os.path.join(savepath, 'FP_grid_distribution.pdf'), transparent=True)
print('threshold:', np.median(grid_values))
print('epsilon:', np.percentile(grid_values, 95))
