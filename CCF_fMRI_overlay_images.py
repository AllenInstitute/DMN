#!/shared/utils.x86_64/python-2.7/bin/python
"""
Generate .png images of CCF with fMRI mask overlaid and white background
"""

import os
import scipy.ndimage
import nrrd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import matplotlib.pyplot as plt
import numpy as np
import platform

if platform.system == 'Linux':
    fmri_datpath = r'DMN paper/fMRI data from Alessandro/fMRI_ICA_map'
    savepath = r'DMN paper/fMRI data from Alessandro/fMRI_ICA_map/fMRI CCF overlap images'
else:
    fmri_datpath = r'DMN paper\fMRI data from Alessandro\fMRI_ICA_map'
    savepath = r'DMN paper\fMRI data from Alessandro\fMRI_ICA_map\fMRI CCF overlap images'

d, _ = nrrd.read(os.path.join(fmri_datpath, 'dmn_mask_z_1_allen_masked_sym.nrrd'))
d17, _ = nrrd.read(os.path.join(fmri_datpath, 'dmn_mask_z_1.7_allen_masked_sym.nrrd'))
combined_map = d + d17

mcc = MouseConnectivityCache(resolution = 25, manifest_file='connectivity\manifest.json')
template, _ = mcc.get_template_volume()
annotation, _ = mcc.get_annotation_volume()
mask = np.where(annotation == 0)
template = template.astype(float)
template[mask] = np.nan

cm = scipy.ndimage.zoom(combined_map, 4, order=0)
cm[mask] = np.nan

for ii in xrange(cm.shape[0]):
    plt.imshow(template[ii], cmap='gray', vmax=255)
    plt.imshow(cm[ii], cmap = 'hot', alpha = 0.5, vmin = 0, vmax = 4)
    plt.axis('off')
    path = os.path.join(savepath, '%s.png' %ii)
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
