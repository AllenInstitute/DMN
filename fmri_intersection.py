#!/shared/utils.x86_64/python-2.7/bin/python

#General imports
from __future__ import division, print_function, absolute_import
import os
import pandas as pd
import sys
import warnings
import nrrd
import numpy as np
        
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

from nileg_utilities import filesystem

#==============================================================================
#
#==============================================================================


class FmriIntersect(object):


    def fmri_mask_path(self, value):
        '''Path to fmri mask
        
        Parameters
        ----------
        value : int
            fmri value
        
        '''
        
        return os.path.join(self.mask_dir, 'fmri_{0}_mask.nrrd'.format(value))
        
    
    def structure_path(self, structure_id):
        '''Path to structure mask
        
        Parameters
        ----------
        structure_id : int
            specifies structure
        
        '''
        
        return os.path.join(
            self.structure_dir, 'structure_{0}_mask'.format(structure_id)
            )
            
    
    def intersection_path(self, structure_id, fmri_value):
        '''Path to fmri-structure intersection 
        mask
        
        Parameters
        ----------
        structure_id : int
            specifies structure
        fmri_value : int
            specifies voxels in fmri data
        
        '''
        
        return os.path.join(
            self.intersection_dir, 
            'structure_{0}_fmri_{1}_mask'.format(structure_id, fmri_value)
            )
        

    def __init__(self, top_level, fmri_path, resolution=100):
        '''Given a ccf-aligned volume with integer labels, 
        generate masks that represent the intersection of those
        labels with ccf structures
        
        Parameters
        ----------
        top_level : str
            directory where results will be stored. Will be
            created if necessary.
        fmri_path : str 
            path to a ccf-shaped nrrd file with integer labels
        resolution : int, optional
            in microns; defaults to 100
        model_dir : str, optional
            this directory should contain a nrrd file called
            'annotation_{resolution}'. This will be used as the 
            ccf annotation. Defaults to the latest model in 
            /projects
        
        '''
    
        self.fmri_path = fmri_path
        
        self.resolution = resolution
        self.top_level = top_level
        filesystem.safe_makedirs(self.top_level)
        self.mask_dir = os.path.join(self.top_level, 'fmri_masks')
        filesystem.safe_makedirs(self.mask_dir)
        self.structure_dir = os.path.join(self.top_level, 'structure_masks')
        filesystem.safe_makedirs(self.structure_dir)
        self.intersection_dir = os.path.join(
            self.top_level, 'structure_intersections'
            )
        filesystem.safe_makedirs(self.intersection_dir)
        self.mcc = MouseConnectivityCache(resolution=self.resolution, 
                                 manifest_file='..\mouse_connecitivity_manifest.json')
        
        
    def make_fmri_masks(self):
        '''Finds unique values in fmri volume and 
        generates a mask for each
        '''
    
        fmri_data = nrrd.read(self.fmri_path)[0]
        
        self.unique_fmri_values = []
        for val in np.unique(fmri_data):
            ival = int(val)
            if abs(ival - val) > sys.float_info.epsilon:
                warnings.warn(
                    'casting {0} to {1} loses information'.format(val, ival), 
                    UserWarning
                    )
            self.unique_fmri_values.append(ival)
            print(self.unique_fmri_values)
        
        for val in self.unique_fmri_values:
            val_mask = np.zeros(fmri_data.shape)
            val_mask[np.where(fmri_data == val)] = 1
            nrrd.write(self.fmri_mask_path(val), val_mask)
            
            
    def get_structure_mask(self, strid, force_update = False):
        '''Load a structure mask
        
        Parameters
        ----------
        structure_ids : list of ints
            specifies the structures
        force_update : bool, optional
            if True -> structure mask will be computed 
            regardless of whether it can be loaded. Defaults 
            to False
            
        Returns
        -------
        ndarray : 
            structure mask. 1 where the structure is present, 0
            elsewhere. from MouseConnectivityCache
        
        '''
        
        if not force_update:
            try:
                return nrrd.read(self.structure_path(strid))[0]
            except (IOError, nrrd.NrrdError) as err:
                try:
                    mask = self.mcc.get_structure_mask(strid)[0]
                    nrrd.write(self.structure_path(strid), mask)
                    return mask
                except (IOError, nrrd.NrrdError) as err:
                    pass
        else:
            mask = self.mcc.get_structure_mask(strid)[0]
            nrrd.write(self.structure_path(strid), mask)
            return mask
        
    def many_structure_masks(self, structure_ids, force_update = False):
        '''Download and save a whole bunch of structure masks
        
        Parameters
        ----------
        structure_ids : list of ints
        
        force_update : Bool
            generate new masks even if new masks exist
        
        '''
        
        for stid in structure_ids:
            self.get_structure_mask(stid, force_update)
            
            
    def get_intersection_mask(self, 
        structure_id, fmri_value, force_update=False
        ):
        '''Loads an intersection mask
        
        Parameters
        ----------
        structure_id : int
            specifies structure
        fmri_value : int
            specifies voxels in fmri data
        force_update : bool, optional
            set to True if you want to compute 
            a new mask regardless of whether one is 
            loadable
            
        '''
        
        if not force_update:
            try:
                return nrrd.read(
                    self.intersection_path(structure_id, fmri_value)
                    )[0]
            except (IOError, nrrd.NrrdError) as err:
                return self.make_intersection_mask(structure_id, fmri_value)
        else:
            return self.make_intersection_mask(structure_id, fmri_value)
            
            
    def make_intersection_mask(self, structure_id, fmri_value):
        '''Compute and store a binary mask that 
        indicates the voxelwise intersection of a 
        structure and fmri value
        
        Parameters
        ----------
        structure_id : int
            specifies structure
        fmri_value : int
            specifies voxels in fmri data
        
        '''
        
        structure_mask = self.get_structure_mask(structure_id)
        fmri_mask, _ = nrrd.read(self.fmri_mask_path(fmri_value))
        
        try:
            intersection = np.multiply(structure_mask, fmri_mask)
            nrrd.write(
            self.intersection_path(structure_id, fmri_value), intersection
            )
            return intersection
        except (TypeError): #This will occur if the structure mask does
            #not exist at the selected resolution
            print(structure_id)
            pass
        
        
    
    def many_intersection_masks(self, structure_ids, fmri_values):
        '''Compute and save multiple fmri-structure 
        intersection masks
        
        Parameters
        ----------
        structure_ids : list of int
            make masks for these structures
        fmri_values : list of int
            make masks for these values
        
        '''
        
        for stid in structure_ids:
            for fmv in fmri_values:
                self.make_intersection_mask(stid, fmv) 

