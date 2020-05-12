# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 18:36:31 2019

@author: jenniferwh
"""
import os
import numpy as np
import pandas as pd
import nrrd
from anatomy.anatomy_api import AnatomyApi
import scipy.ndimage as ndimage

path = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\Figure 6 target-defined matrix'
maskpath = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\fMRI masks'
savepath = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\Data files'

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
resolution = 100 #100 um for fMRI overlap
mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                            resolution=resolution)
# Downsampled CAV grids are stored in my network folder
cav_basepath = r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN_paper'
cav_dir = os.path.join(cav_basepath, 'downsampled_CAV', 'grid') # Note the rabies polylgons are in this directory as well
structure_dir = os.path.join(cav_basepath, 'structure_masks')

def make_fmri_masks(fmri_mask, savepath):
    """Finds unique values in fmri volume and generates a mask for each.
    """
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    unique_fmri_values = []
    for val in np.unique(fmri_mask):
        ival = int(val)
        unique_fmri_values.append(ival)
    
    for val in unique_fmri_values:
        val_mask = np.zeros(fmri_mask.shape)
        val_mask[np.where(fmri_mask == val)] = 1
        nrrd.write(os.path.join(
                savepath, 
                'fmri_{0}_mask.nrrd'.format(val)), val_mask)

dmn_mask, _ = nrrd.read(os.path.join(maskpath, 'dmn_mask_and_core.nrrd'))
mask_dir = os.path.join(path, 'fmri_masks')
make_fmri_masks(dmn_mask, mask_dir)
fmri0, _ = nrrd.read(os.path.join(mask_dir, 'fmri_0_mask.nrrd'))
fmri2, _ = nrrd.read(os.path.join(mask_dir, 'fmri_2_mask.nrrd'))

aapi = AnatomyApi()
CAVis = aapi.get_image_series_by_workflow([471789262, 304950893])
td_dataset = pd.read_csv(r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\target_defined_dataset.csv')

good_experiments = td_dataset[td_dataset['include'] != 'no']['image_series_id'].unique()
CAVis = [isid for isid in CAVis if isid in good_experiments]
print(len(CAVis))

for isid in CAVis:
    if not os.path.isdir(os.path.join(cav_dir, str(isid))):
        limspath = os.path.join(aapi.get_storage_directory(isid), 'grid')
        gridfile, _ = nrrd.read(os.path.join(limspath, 'cav_density_10.nrrd'))
        dsgrida = ndimage.zoom(gridfile, 0.4)
        newpath = os.path.join(cav_dir, str(isid))
        if not os.path.isdir(newpath):
            os.mkdir(newpath)
        nrrd.write(os.path.join(newpath, 'CAV_density_25.nrrd'), dsgrida)   
        dsgridb = ndimage.zoom(gridfile, 0.1)
        nrrd.write(os.path.join(newpath, 'CAV_density_100.nrrd'), dsgridb)       
        
polysize = []
OL0 = []
OL2 = []
CAVproj_path = os.path.join(cav_basepath, 'downsampled_CAV', 'segmentation_in_CAV_mask')
for isid in CAVis:
    cavpath = os.path.join(cav_dir, str(isid))
    CAV, _ = nrrd.read(os.path.join(cavpath, 'cav_density_100.nrrd'))
    gridpath = os.path.join(aapi.get_storage_directory(isid), 'grid')
    proj, _ = nrrd.read(os.path.join(gridpath, 'projection_density_100.nrrd'))
    CAVmask = np.zeros(np.shape(CAV))
    CAVmask[np.where(CAV>0.001)] = 1
    CAVproj = proj * CAVmask
    projpath = os.path.join(CAVproj_path, str(isid))
    if not os.path.exists(projpath):
        os.mkdir(projpath)
    nrrd.write(os.path.join(projpath, 'projections_in_CAV_100.nrrd'), CAVproj)
    polysize.append(np.sum(CAVproj))
    OL0.append(np.sum(fmri0 * CAVproj))
    OL2.append(np.sum(fmri2 * CAVproj))

CAV_dat = pd.DataFrame({'image_series_id': CAVis, 'CAV_size': polysize, 'out_DMN': OL0, 'in_DMN': OL2})
CAV_dat['CAV_percent_DMN'] = CAV_dat[['in_DMN']].sum(axis = 1)/CAV_dat['CAV_size']*100

# %% Find target injection structures
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
resolution = 25 # use 25 um resolution for injection structures
mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                            resolution=resolution)
summary_structure_ids = aapi.get_summary_structure_data('id')

def get_structure_masks(structure_id, savepath=None):
    '''Get structure masks from Allenllen API
    
    Parameters
    ----------
    structure_id : int
        specifies structure
    savepath : string
        location to store masks for futre use
    
    '''
    structure_mask, _ = mcc.get_structure_mask(structure_id)
    if savepath:
        if not os.path.exists(savepath):
            os.mkdir(savepath)
            nrrd.write(os.path.join(savepath, 'structure_{0}_mask.nrrd'.format(structure_id)), structure_mask)
    else:
        return structure_mask

structure_dir_resolution = os.path.join(structure_dir, str(resolution))
for structure_id in summary_structure_ids:
    get_structure_masks(structure_id, structure_dir_resolution)
    
CAVproj_path = os.path.join(cav_basepath, 'downsampled_CAV', 'segmentation_in_CAV_mask')
for isid in CAVis:
    cavpath = os.path.join(cav_dir, str(isid))
    CAV, _ = nrrd.read(os.path.join(cavpath, 'cav_density_25.nrrd'))
    gridpath = os.path.join(aapi.get_storage_directory(isid), 'grid')
    proj, _ = nrrd.read(os.path.join(gridpath, 'projection_density_25.nrrd'))
    CAVmask = np.zeros(np.shape(CAV))
    CAVmask[np.where(CAV>0.001)] = 1
    CAVproj = proj * CAVmask
    projpath = os.path.join(CAVproj_path, str(isid))
    if not os.path.exists(projpath):
        os.mkdir(projpath)
    nrrd.write(os.path.join(projpath, 'projections_in_CAV_25.nrrd'), CAVproj)

# %% Get primary and secondary injection structures for CAV data. 
# TODO: Find out how to get CAV grids through API/SDK and update
for isid in CAVis:
    CAV, _ = nrrd.read(os.path.join(CAVproj_path, str(isid), 'projections_in_CAV_25.nrrd'))
    maxOL = 0
    OL = 0
    secondary = 997
    primary = 997
    for mask in os.listdir(os.path.join(structure_dir, '25')):
        structure = mask[10:-10]
        if not any(structure in parent_structures for parent_structures in ['997', '688', '315', '8', '1080', '1009', '998']): 
            #parent structures: root, cortex, isocortex, basic, hippocampus, fiber tracts 
            strmask, _ = nrrd.read(os.path.join(structure_dir, '25', mask))
            OL = np.sum(strmask * CAV)
            if OL > maxOL:
                maxOL = OL
                secondary = primary
                primary = mask[10:-10]
    CAV_dat.loc[CAV_dat['image_series_id'] == isid, 'primary target'] = int(primary)
    CAV_dat.loc[CAV_dat['image_series_id'] == isid, 'secondary target'] = int(secondary)

structure_tree = mcc.get_structure_tree()
ai_map = structure_tree.get_id_acronym_map()
ia_map = {value:key for key, value in ai_map.items()}
CAV_dat['primary target acronym'] = [ia_map[primary] for primary in CAV_dat['primary target']]
CAV_dat['secondary target acronym'] = [ia_map[secondary] for secondary in CAV_dat['secondary target']]

# %% Add source injection structure information
td_experiments = pd.DataFrame(mcc.get_experiments(cre=['Ai75(RCL-nt)']))
td_experiments.rename(columns = {'id': 'image_series_id'}, inplace = True)
td_experiments = CAV_dat.merge(td_experiments[['gender', 'image_series_id', 'injection_structures', 'injection_volume', 
                                               'injection_x', 'injection_y', 'injection_z',
                                               'primary_injection_structure', 'structure_abbrev',
                                               'structure_name']], on='image_series_id', how='left')
#%% GLM for DMN coefficient
def shrink_to_mask(x, mask):
    """
    returns x in space of the mask
    
    NOTE: flattened!
    
    Parameters
    ----------
    x : array-like
        array to shrink (sample)
        could be array or mask itself
    mask : array-like, SAME SHAPE AS X
        mask with which to shrink (sample)
        
    Returns
    -------
    flattened elements of x in mask
    """
    if not all(np.equal(x.shape, mask.shape)):
        raise ValueError("x and mask must have same shape!")
    
    return x[np.where(mask)]

from scipy.spatial.distance import cdist

def compute_distances_from_centroid(centroid, mask):
    '''computes euclidean distance from centroid to all points in given mask
    Parameters
    ----------
      centroid : list or np.array
        3d coordinates of centroid in units of voxels! May have to scale from 
        output of mca.calculate_injection_centroid() by dividing by resolution
      mask : np.array
        binary mask of region in which the distances are wanted
    
    Returns
    -------
      distances : list
        list of distances from every voxel in masked region to centroid
    '''
    
    mask_coordinates = np.argwhere(mask)
    centroid_ = np.atleast_2d(centroid)
    
    distances = cdist(centroid_, mask_coordinates, metric='euclidean')
    return distances.ravel()

from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
mca = MouseConnectivityApi()

def calculate_injection_centroid(injection_density,
                                     injection_fraction,
                                     resolution=1):
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

  #%%

import statsmodels.api as sm  
def fit_glm_CAV_removed(categorical_var, distances, projections):
    '''inputs:
    1. categorical variable (in_or_out)
    2. distances
    3. projections
    #(x1 is distance coeff, x2 is "dmn" coeff)
    # fit glm for each experiment
    '''
    coeff1 = [] #distance
    coeff2 = [] #DMN
    tvals = []
    pvals = []
    for exp in range(len(distances)):
        groups = np.array(categorical_var[exp])

        dummy = sm.categorical(groups, drop=True)
        x = distances[exp]

        # drop reference category
        X = np.column_stack((x, dummy[:,1:]))
        X = sm.add_constant(X, prepend=False)

        # y Use log projection density
        y = np.log10 ( projections[exp] + 0.0000001 )
    
        # fit
        fit = sm.OLS(y, X).fit()
    
        # add coeff
        coeff1 += [fit.params[0]]
        coeff2 += [fit.params[1]]
        tvals += [fit.tvalues]
        pvals += [fit.pvalues]
    return coeff1, coeff2, tvals, pvals  

#%%
cav_basepath = r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN_paper'
# calculate centroids & distances
centroids = []
ratios_in = []
projections_in = []
distances = []
projections = []
in_or_out_masks = []

resolution = 100
mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                            resolution=resolution)
iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
iso_mask, _ = mcc.get_structure_mask(iso['id']) # 100 um resolution

# calculate centroids & distances
for exp in td_experiments['image_series_id']:
    try:
        inj_frac = mcc.get_injection_fraction(exp)[0]
        inj_den = mcc.get_injection_density(exp)[0]
        proj_den = mcc.get_projection_density(exp)[0]
        data_mask = mcc.get_data_mask(exp)[0]
    except:
        grid_path = os.path.join(aapi.get_storage_directory(exp), 'grid')
        inj_frac, _ = nrrd.read(os.path.join(grid_path, 'injection_fraction_100.nrrd'))
        inj_den, _ = nrrd.read(os.path.join(grid_path, 'injection_density_100.nrrd'))
        proj_den, _ = nrrd.read(os.path.join(grid_path, 'projection_density_100.nrrd'))
        data_mask, _ = nrrd.read(os.path.join(grid_path, 'data_mask_100.nrrd'))
    
    CAV_dir = os.path.join(cav_basepath, 'downsampled_CAV', 'grid', str(exp))
    cav_den, _ = nrrd.read(os.path.join(CAV_dir, 'CAV_density_100.nrrd'))
    
    # extend isocortex mask to exclude CAV injection site and data mask. 
    # This mask will be used to index all data
    iso_cav_removed = np.copy(iso_mask)
    iso_cav_removed[np.where(cav_den > 0.001)] = 0
    iso_cav_removed[np.where(data_mask < 0.5)] = 0

    # create unique in_jor_out mask per experiment
    in_or_out = shrink_to_mask(dmn_mask, iso_cav_removed)
    
    # apply data mask. Projection density was already masked when the CAV signal was removed with 
    # remove_CAV_projections.py
    inj_den = np.multiply(inj_den, data_mask)
    inj_den[np.where(cav_den>0.001)] = 0 #This should not happen except for overlapping injections
    inj_den[np.where(data_mask < 0.5)] = 0
    inj_frac = np.multiply(inj_frac, data_mask)
    
    # compute normalized projection density
    proj_den_cav_removed = np.copy(proj_den)
    proj_den_cav_removed[np.where(cav_den > 0.001)] = 0
    proj_den_cav_removed[np.where(data_mask < 0.5)] = 0
    npd = proj_den_cav_removed / inj_den.sum()

    # centroid (res = 1 so as to keep in voxel space)
    centroid = mca.calculate_injection_centroid(inj_den,
                                                inj_frac,
                                                resolution = 1)
    
    # compute distances
    distance = compute_distances_from_centroid(centroid,iso_cav_removed)
    
    # find ratio inside mask
    ratio = sum(inj_den[np.where(dmn_mask)]) / sum(inj_den.flatten())
    projection_ratio = sum(proj_den_cav_removed[np.where(dmn_mask)]) / sum(proj_den_cav_removed.flatten())

    # add to lists
    centroids += [centroid]
    ratios_in += [ratio]
    projections_in += [projection_ratio]
    distances += [distance]
    projections += [ npd[np.where(iso_cav_removed)] ]
    in_or_out_masks.append( in_or_out )
    
d_coeff, dmn_coeff, tvals, pvals  = fit_glm_CAV_removed(in_or_out_masks, distances, projections)

td_experiments['injection dmn fraction CAV removed'] = ratios_in
td_experiments['projection dmn fraction CAV removed'] = projections_in
td_experiments['distance coefficient CAV removed'] = d_coeff
td_experiments['DMN coefficient CAV removed'] = dmn_coeff
td_experiments['t values CAV removed'] = tvals
td_experiments['p values CAV removed'] = pvals

#%%
td_experiments.to_csv(os.path.join(
        savepath, 
        'CAV_DMN_data.csv'), index = False)
