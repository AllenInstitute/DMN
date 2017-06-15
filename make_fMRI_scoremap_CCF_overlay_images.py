# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:39:52 2017

@author: jenniferwh
"""

import os
import scipy.ndimage
import nrrd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import matplotlib.pyplot as plt

#turn interactive plotting off
plt.ioff()

fmri_datpath = r'\\AIBSDATA\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\fMRI data from Alessandro\fMRI_ICA_map'
savepath = r'\\AIBSDATA\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\fMRI data from Alessandro\fMRI_ICA_map\fMRI CCF overlap images'
if not os.path.exists(savepath):
    os.mkdir(savepath)
  
#load two fMRI masks: core regions (z score > 1.7) and standard (z score > 1 (d))  
d, _ = nrrd.read(os.path.join(fmri_datpath, 'dmn_mask_z_1_allen_masked_sym.nrrd'))
d17, _ = nrrd.read(os.path.join(fmri_datpath, 'dmn_mask_z_1.7_allen_masked_sym.nrrd'))

mcc = MouseConnectivityCache(resolution = 25, manifest_file='connectivity/manifest.json')
template, template_info = mcc.get_template_volume()

#Resample fMRI masks at 25 um
d25 = scipy.ndimage.zoom(d, 4, order=0)
d1725 = scipy.ndimage.zoom(d17, 4, order=0)

#Combine two masks into one for plotting. Had to mess with these numbers to get colors I liked
masks = d25
masks[d25 > 0] = 25
masks[d1725 > 0] = 75
d1725[d1725 > 0] = 1

fig, ax = plt.subplots(1,1, figsize=(10, 10))
for ii in xrange(d25.shape[0]):
    slc = d25[ii]
    ax.imshow(template[ii], cmap='gray', vmin=template.min(), vmax=template.max())
    ax.imshow(d1725[ii], cmap='Purples', alpha=0.3)
    ax.imshow(masks[ii], cmap = 'nipy_spectral', alpha = 0.3)
    plt.axis('off')
    path = os.path.join(savepath, '%s.png' %ii)
    plt.savefig(path, dpi=200, bbox_inches='tight')