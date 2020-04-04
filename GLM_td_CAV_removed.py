# -*- coding: utf-8 -*-
"""
Created on Wed Aug 02 16:18:34 2017

@author: jenniferwh
"""
import nrrd
import os
import numpy as np
import pandas as pd
import platform
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist
import statsmodels.api as sm
from anatomy.anatomy_api import AnatomyApi
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
aapi = AnatomyApi()
mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                            resolution=100)
mca = MouseConnectivityApi()
structure_tree = mcc.get_structure_tree()

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
    
def compute_distance_from_centroid(centroid, mask):
    
    # get coordinates for all (returns [ [a,b,..,z],[a,b,...,z],[a,b,...,z] ] format)
    coordinates = np.where( mask )

    # but we want in [ [a,b,c],[a,b,c],...,[a,b,c] ] format
    coordinates = np.stack(coordinates, axis = -1)
    
    # compute distances between centroid and all other points in mask
    # (returns array(array(distances)))
    distances = cdist(np.atleast_2d(centroid),coordinates)

    return distances.ravel()

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
        y = np.log10 ( projections[exp] + 1 )
    
        # fit
        fit = sm.OLS(y, X).fit()
    
        # add coeff
        coeff1 += [fit.params[0]]
        coeff2 += [fit.params[1]]
        tvals += [fit.tvalues]
        pvals += [fit.pvalues]
    return coeff1, coeff2, tvals, pvals
    
# want only inside isocortex
iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
iso_mask = mcc.get_structure_mask(iso['id'])[0]

# grab some experiments
td_experiments = pd.read_csv(r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\target_defined_dataset.csv')
td_experiments = td_experiments[td_experiments['include'] == 'yes']
print(len(td_experiments))

# DMN maksks
if platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN'
elif platform.system() == 'Darwin':
    path = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN'
outpath = os.path.join(path, 'data_files')
masks, _ = nrrd.read(os.path.join(path, 'fMRI_masks', 'dmn_mask_and_core.nrrd'))
dmn_mask = np.zeros(masks.shape)
dmn_mask[np.where(masks > 0)] = 1
core_mask = np.zeros(masks.shape)
core_mask[np.where(masks == 2)] = 1

# calculate centroids & distances
centroids = []
inj_ratios = []
proj_ratios = []
core_inj_ratios = []
core_proj_ratios = []
distances = []
projections = []
in_or_out_masks = []
in_or_out_core_masks = []
cav_basepath = r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN_paper'
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
    
    inj_den = np.multiply(inj_den, data_mask)
    CAV_dir = os.path.join(cav_basepath, 'downsampled_CAV', 'grid', str(exp))
    cav_den, _ = nrrd.read(os.path.join(CAV_dir, 'CAV_density_100.nrrd'))
    
    # extend isocortex mask to exclude CAV injection site and data mask. 
    # This mask will be used to index all data
    iso_cav_removed = np.copy(iso_mask)
    iso_cav_removed[np.where(cav_den > 0.001)] = 0
    iso_cav_removed[np.where(data_mask < 0.5)] = 0

    # create unique in_or_out mask per experiment
    in_or_out_CAV = shrink_to_mask(dmn_mask, iso_cav_removed)
    in_or_out_CAV_core = shrink_to_mask(core_mask, iso_cav_removed)
    
    # apply data mask. Projection density was already masked when the CAV signal was removed with 
    # remove_CAV_projections.py
    inj_den = np.multiply(inj_den, data_mask)
    inj_den[np.where(cav_den>0.001)] = 0 #This should not happen except for overlapping injections
    inj_den[np.where(data_mask < 0.5)] = 0

    # centroid (res = 1 so as to keep in voxel space)
    centroid = mca.calculate_injection_centroid(inj_den,
                                                inj_frac,
                                                resolution = 1)
    
    # compute distances
    distance = compute_distance_from_centroid(centroid,iso_cav_removed)
    
    proj_den_cav_removed = np.copy(proj_den)
    proj_den_cav_removed[np.where(cav_den > 0.001)] = 0
    proj_den_cav_removed[np.where(data_mask < 0.5)] = 0
    proj_iso = proj_den_cav_removed * iso_mask
    
    # find ratio inside mask
    inj_ratio = sum(inj_den[np.where(dmn_mask)]) / sum(inj_den.flatten())
    core_inj_ratio = sum(inj_den[np.where(core_mask)]) / sum(inj_den.flatten())
    proj_ratio = sum(proj_iso[np.where(dmn_mask)]) / sum(proj_iso.flatten())
    core_proj_ratio = sum(proj_iso[np.where(core_mask)]) / sum(proj_iso.flatten())

    # add to lists
    centroids += [centroid]
    inj_ratios += [inj_ratio]
    core_inj_ratios += [core_inj_ratio]
    proj_ratios += [proj_ratio]
    core_proj_ratios += [core_proj_ratio]
    distances += [distance]
    projections += [ proj_den[np.where(iso_cav_removed)] ]
    in_or_out_masks.append( in_or_out_CAV )
    in_or_out_core_masks.append( in_or_out_CAV_core )

# fit glm
d_coeff, dmn_coeff, tvals, pvals  = fit_glm_CAV_removed(
    in_or_out_masks, distances, projections)
d_coeff_core, dmn_coeff_core, tvals_core, pvals_core  = fit_glm_CAV_removed(
    in_or_out_core_masks, distances, projections)

x=inj_ratios
y=dmn_coeff
dmn_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(dmn_fit.summary())

x=inj_ratios
y=dmn_coeff_core
core_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(core_fit.summary())

x=inj_ratios
y=d_coeff
distance_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(distance_fit.summary())

ctx_glm_dat = pd.DataFrame({'id': td_experiments['image_series_id'], 
                       'injection structure': td_experiments['source'],
                       'injection dmn fraction': inj_ratios,
                       'projection dmn fraction': proj_ratios,
                       'distance coefficient': d_coeff,
                       'DMN coefficient': dmn_coeff,
                       'DMN t values': tvals,
                       'DMN p values': pvals,
                       'injection core fraction': core_inj_ratios,
                       'projection core fraction': core_proj_ratios,
                       'core distance coefficient': d_coeff_core,
                       'DMN core coefficient': dmn_coeff_core,
                       'core t values': tvals_core,
                       'core p values': pvals_core})
ctx_glm_dat.to_csv(os.path.join(outpath, 'td_ctx_injections_DMN_and_core_projections_coefficients.csv'),
              index = False)

x=inj_ratios
y=dmn_coeff
dmn_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(dmn_fit.summary())
