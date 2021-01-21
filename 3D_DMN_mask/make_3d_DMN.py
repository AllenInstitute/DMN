#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 19:19:17 2017

@author: jen
"""
import sys
from mayavi import mlab
from tvtk.api import tvtk
from tvtk.util.ctf import ColorTransferFunction
from scipy.ndimage.filters import gaussian_filter
import nrrd
import os
import platform
import scipy.ndimage

# get a binary mask
if platform.system() == 'Windows':
    basepath = r'fMRI data from Alessandro\fMRI_ICA_map'
if platform.system() == 'Darwin':
    basepath = r'fMRI data from Alessandro/fMRI_ICA_map'

mask1, _ = nrrd.read(os.path.join(basepath, 'dmn_mask_z_1_allen_masked_sym.nrrd'))
mask1_25um = scipy.ndimage.zoom(mask1, 4, order=0)
nrrd.write(os.path.join(basepath, 'dmn_mask_z_1_allen_masked_sym_25um.nrrd'), mask1_25um, options={ 'encoding': 'raw' })
mask2, _ = nrrd.read(os.path.join(basepath, 'dmn_mask_z_1.7_allen_masked_sym.nrrd'))
mask2_25um = scipy.ndimage.zoom(mask2, 4, order=0)
nrrd.write(os.path.join(basepath, 'dmn_mask_z_1.7_allen_masked_sym_25um.nrrd'), mask2_25um, options={ 'encoding': 'raw' })

template, _ = nrrd.read(os.path.join(basepath, 'template_white_background.nrrd'))

gaussian_mask1 = gaussian_filter(mask1_25um.astype(float), sigma=3)
gaussian_mask2 = gaussian_filter(mask2_25um.astype(float), sigma=3)
gaussian_template = gaussian_filter(template.astype(float), sigma = 0.3)
template_field = mlab.pipeline.scalar_field(gaussian_template)
mask1_field = mlab.pipeline.scalar_field(gaussian_mask1)
mask2_field = mlab.pipeline.scalar_field(gaussian_mask2)

ctf = ColorTransferFunction()
ctf.add_rgb_point(0.5, 1, 1, 1)

template = mlab.pipeline.volume(template_field)
DMN = mlab.pipeline.volume(mask2_field)
DMN._volume_property.set_color(ctf)
DMN._ctf = ctf
DMN.update_ctf = True
DMNcore = mlab.pipeline.volume(mask1_field)
mlab.show()
