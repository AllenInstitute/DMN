#!/shared/utils.x86_64/python-2.7/bin/python
from __future__ import division, print_function, absolute_import
import os
import json
import argparse

from allensdk.internal.core.lims_utilities import safe_system_path

from nileg_projects.fmri_unionization.fmri_intersection import FmriIntersect
from nileg_projects.fmri_unionization.fmri_signal_unionization import \
    FmriUnionize, FmriUnionizeBlock
    
from nileg_utilities.lims_queries import \
    passed_from_image_series_id, image_series_ids_from_project_code, \
    image_series_ids_from_workflow_name
from nileg_utilities.filesystem import safe_makedirs
from nileg_utilities.ontology import StructureArray

'''
Create two folders in the current directory called "inputs" and "outputs".

Inside the "inputs" folder there should be two folders, one called "configs" 
and one called "scoremaps". Inside the configs folder should be a json file
that specifies the structures, image series ids, fMRI mask, and resolution to
use. Inside the "scoremaps" folder are the fMRI masks in nrrd format.

If the structure masks, structure interesections, and fMRI masks have already
been created, they should be placed in folders in the "outputs" folder in 
three seperate folders named "structure_masks", "structure_intersections",
and "fmri_masks". If these have not been created, they will be made with this
script. Inside the "outputs" folder, a folder will be created for each dataset.
This filename is specified by the input to "top level"

Inputs:
1. config_path. This is the name of the json file inside the configs folder
that should be used (inputs\configs will be added to path)
2. top_level. This is the name of the outputs folder that will be created
for this dataset.
3. Make structure intersections. These structure x fMRI overlap masks can
be reused for different datasets. Default = False to use existing structure
intersections, set to True to make new structure x fMRI mask intersections

'''
def main():
    
    basepath = os.getcwd()
    config_path = safe_system_path(os.path.join(
        basepath, 'inputs', 'configs', args.config_path))
    with open(config_path, 'rb') as cfgf:
        config = json.load(cfgf)
    
    fmri_values = [0, 1, 2]
    base_sm_path = os.path.join(basepath, 'inputs')
    base_top_level = os.path.join(basepath, 'outputs')
    
    top_level = safe_system_path(os.path.join(base_top_level, args.top_level))
    safe_makedirs(top_level)
    sm_path = safe_system_path(os.path.join(base_sm_path, 'scoremaps', 
                                            config['fmri_mask']))
    if args.make_structure_intersections:
        inter = FmriIntersect(
            top_level=base_top_level, 
            fmri_path=sm_path, 
            resolution=config['resolution'],
            model_dir = r'E:\fMRI_unionize\annotation\ccf_2016'
            )
        inter.make_fmri_masks()
        inter.many_structure_masks(*config['structure_ids'])
        inter.many_intersection_masks(config['structure_ids'], fmri_values)
    
        mask_dir = os.path.join(base_top_level, 'structure_intersections')
    else:
        mask_dir = safe_system_path(os.path.join(base_top_level,
                                                  'structure_intersections'))
        
    fbun = FmriUnionizeBlock(
        top_level= os.path.join(top_level, 'block'),
        mask_dir = mask_dir, 
        plaque_dataset=bool(config['plaque_dataset']),
        resolution=config['resolution'],
        store=False
        )  

    fbun.unionize_many_series( 
        config['image_series_ids'], fmri_values, config['structure_ids']
        )
        

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('config_path', type=str)
    parser.add_argument('top_level', type=str)
    parser.add_argument('make_structure_intersections', type=bool, default=False)
    
    args = parser.parse_args()

    main()
