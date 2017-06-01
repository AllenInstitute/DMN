# -*- coding: utf-8 -*-
"""
Created on Fri May 26 16:45:44 2017

@author: jenniferwh
"""

import numpy as np
import os
import nrrd
from anatomy.anatomy_api import AnatomyApi
import pandas as pd

savepath = r'\\AIBSDATA\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN paper\consensus targets'
source = 'TEa_ECT_ENTl'

def calculate_injection_centroid(injection_density,
                                     injection_fraction,
                                     resolution=100):
        '''
        Compute the centroid of an injection site.
        
        Parameters
        ----------
        
        injection_density: np.ndarray
            The injection density volume of an experiment

        injection_fraction: np.ndarray
            The injection fraction volume of an experiment

        '''

        # find all voxels with injection_fraction > 0
        injection_voxels = np.nonzero(injection_fraction)
        injection_density_computed = np.multiply(injection_density[injection_voxels],
                                                 injection_fraction[injection_voxels]) 
        sum_density = np.sum(injection_density_computed)
    
        # compute centroid in CCF coordinates
        if sum_density > 0 :
            centroid = np.dot(injection_density_computed,
                              zip(*injection_voxels)) / sum_density * resolution
        else:
            centroid = None
        
        return centroid
def get_gridpaths(isids):
    filepaths = AnatomyApi()
    grid_paths = []
    for imser in isids:
        grid_paths.append(os.path.join(filepaths.get_storage_directory(imser), 'grid'))
    return grid_paths

def get_centroids(grid_paths):
    centroids = []
    for path in grid_paths:
        #Read in injection_density and injection_fraction files
        INJDENS, metaINJDENS = nrrd.read(os.path.join(path, 'injection_density_25.nrrd'));
        INJFRAC, metaINJFRAC = nrrd.read(os.path.join(path, 'injection_fraction_25.nrrd'));
        #Calculate centroid
        centroids.append(calculate_injection_centroid(INJDENS, INJFRAC, 25))
    return centroids

### Keep track of wt data separately because these centroids need to be flipped to the other hemisphere  
wtLIMSIDs = [127397469,146553971,180403712,180435652]
LIMSIDs = [484612961,485729987,525413115,528511967,529428776]
allIDs = wtLIMSIDs + LIMSIDs

wt_paths = get_gridpaths(wtLIMSIDs)
grid_paths = get_gridpaths(LIMSIDs)
centroids = get_centroids(grid_paths)
wt_centroids = get_centroids(wt_paths)

### flip wt data to left hemisphere
wt_centroids_flipped = []
for centroid in wt_centroids:
    flipped_ml = 5700 - (centroid[2] - 5700)
    wt_centroids_flipped.append(np.array([centroid[0], centroid[1], flipped_ml]))
    
all_centroids = wt_centroids_flipped + centroids

df = pd.DataFrame(index = allIDs, columns = allIDs)

for counterA, centroidA in enumerate(all_centroids):
    for counterB, centroidB in enumerate(all_centroids):
        df.ix[allIDs[counterA], allIDs[counterB]] = np.linalg.norm(centroidA - centroidB)
        
df.to_csv(os.path.join(savepath, '%s_centroid_table.csv' %source))