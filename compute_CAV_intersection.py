# -*- coding: utf-8 -*-
"""
Created on Tue Feb 07 17:33:04 2017

@author: jenniferwh
"""

from __future__ import division, print_function, absolute_import
import os
import json
import argparse

from allensdk.internal.core.lims_utilities import safe_system_path

from nileg_projects.fmri_unionization.fmri_intersection import FmriIntersect
from nileg_projects.fmri_unionization.fmri_signal_unionization import \
    FmriUnionize, FmriUnionizeBlock
    
#from nileg_utilities.lims_queries import \
#    passed_from_image_series_id, image_series_ids_from_project_code, \
#    image_series_ids_from_workflow_name
from nileg_utilities.filesystem import safe_makedirs
#from nileg_utilities.ontology import StructureArray

def main():
    
    config_path = safe_system_path(os.path.join(    
        r'\\AIBSDATA\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\fMRI data from Alessandro\fMRI CCF overlap improved registration\CAV unionize\inputs\configs', args.config_path))
    with open(config_path, 'rb') as cfgf:
        config = json.load(cfgf)
    
    fmri_values = [0, 1, 2, 3, 4]
    base_sm_path = r'\\AIBSDATA\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\fMRI data from Alessandro\fMRI CCF overlap improved registration\CAV unionize\inputs\scoremaps'
    base_top_level = r'\\AIBSDATA\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\fMRI data from Alessandro\fMRI CCF overlap improved registration\CAV unionize\outputs'

    top_level = base_top_level
    safe_makedirs(top_level)
    
    sm_path = safe_system_path(os.path.join(base_sm_path, config['fmri_mask']))

    
    inter = FmriIntersect(top_level=top_level, fmri_path=sm_path, resolution=config['resolution'])
    inter.make_fmri_masks()
    inter.many_structure_masks(*config['structure_ids'])
    inter.many_intersection_masks(config['structure_ids'], fmri_values)
    
    mask_dir = os.path.join(top_level, 'structure_intersections')

    fbun = FmriUnionizeBlock(
        top_level= os.path.join(top_level, 'block'),
        mask_dir = mask_dir, 
        resolution=config['resolution'], 
        store=False
        )
    fbun.unionize_many_series( 
        config['image_series_ids'], fmri_values, config['structure_ids']
        )
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('config_path', type=str)
    #parser.add_argument('top_level', type=str)
    #parser.add_argument('mask_dir', type=str, default=False)
    
    args = parser.parse_args()

    main()
