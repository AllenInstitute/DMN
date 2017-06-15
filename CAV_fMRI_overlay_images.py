# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 18:15:19 2017

@author: jenniferwh
"""

import os
import scipy.ndimage as ndimage
import nrrd
import lims_utilities as lu
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import matplotlib.pyplot as plt
import numpy as np
from anatomy.anatomy_api import AnatomyApi

basepath = r'\\AIBSDATA\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\fMRI data from Alessandro\fMRI CCF overlap improved registration'
outpath = os.path.join('fMRI ccf overlap images', 'CAV injection sites')
savepath = os.path.join(basepath, outpath)
if not os.path.exists(savepath):
    os.mkdir(savepath)
    
mcc = MouseConnectivityCache(resolution = 25, manifest_file='connectivity/manifest.json')
template, template_info = mcc.get_template_volume()

d, meta = nrrd.read(os.path.join(basepath, 'scoremap_MRI_allen_ccfv3_25um_masked.nrrd'))

image_series_query = '''
    select im.id from image_series im
    join specimens_workflows sw on sw.specimen_id = im.specimen_id
    where sw.workflow_id = 471789262
    and im.workflow_state like 'passed'
'''
image_series_results = lu.query(image_series_query)
isr = [iser['id'] for iser in image_series_results]

aapi = AnatomyApi()
grid_paths = []
for isid in isr:
    grid_path = os.path.join(aapi.get_storage_directory(isid), 'grid')
    grid_paths.append([isid, grid_path])

print(len(isr))
print(len(grid_paths))
        
for line in grid_paths:
    print line[0]
    ind, ind_info = nrrd.read(os.path.join(line[1], 'cav_density_10.nrrd'))
    ind = ndimage.zoom(ind, 0.4)
    inj_site = []
    for ix in np.unique(np.nonzero(ind)):
        inj_site.append([ix, np.count_nonzero(ind[ix])])
    if inj_site != []:
        maxi = []
        maxs = 0
        for ij in inj_site:
            if ij[1] > maxs:
                maxs = ij[1]
                maxi = ij[0]
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.imshow(template[maxi], cmap='gray', vmin=template.min(), vmax=template.max())
        ax.imshow(d[maxi], cmap='magma', alpha = 0.3)
        ax.imshow(ind[maxi], cmap = 'hot', alpha = 0.3)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        filename = os.path.join(savepath, str(line[0])+'.png')
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()