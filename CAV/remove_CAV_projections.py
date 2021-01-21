# -*- coding: utf-8 -*-
"""
Created on Sun Nov 05 14:18:41 2017

@author: jenniferwh
"""
import os
import nrrd
import numpy as np

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

print('first test')
output_directory = r'DMN paper/remove_CAV_signal/grid'
datpaths = {561511939: u'image_series_561511939/'} #use anatomy_api to get paths instead

mcc = MouseConnectivityCache(resolution = 10)
print('mcc')
structure_tree = mcc.get_structure_tree()
iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
iso_mask = mcc.get_structure_mask(iso['id'])[0]
print('iso_mask')

def remove_CAV_projections():
    
    for imser in [475616128]:
        print(imser)
        grid_path = os.path.join(datpaths[imser], 'grid')
        proj_density, _ = nrrd.read(os.path.join(grid_path, 'projection_density_10.nrrd'))
        CAV_density, _ = nrrd.read(os.path.join(grid_path, 'CAV_density_10.nrrd'))  
        data_mask, _ = nrrd.read(os.path.join(grid_path, 'data_mask_10.nrrd'))
        masked_density = np.multiply(proj_density, data_mask)
        masked_density[np.where(CAV_density > 0)] = np.nan
        iso_masked = np.copy(iso_mask)
        iso_masked[np.where(CAV_density > 0)] = np.nan
        out_dir = os.path.join(
            output_directory, '{0}'.format(imser)
            )
        print(out_dir)
        output_mask_path = os.path.join(out_dir, 'iso_masked_10.nrrd')
        output_nrrd_path = os.path.join(out_dir, 'masked_density_10.nrrd')
        nrrd.write(output_mask_path, iso_masked)
        nrrd.write(output_nrrd_path, masked_density)

remove_CAV_projections()
