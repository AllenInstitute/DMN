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
from scipy.spatial.distance import cdist
import statsmodels.api as sm
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi

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

def calculate_centroids_and_distances(experiments, functional_mask, structural_mask):
    # calculate centroids & distances
    centroids = []
    ratios_in = []
    distances = []
    projections = []
    projection_ratios = []

    for exp in experiments['id']:
        print(exp)
        # injection/projection data
        inj_frac = mcc.get_injection_fraction(exp)[0]
        inj_den = mcc.get_injection_density(exp)[0]
        proj_den = mcc.get_projection_density(exp)[0]
        data_mask = mcc.get_data_mask(exp)[0]

        inj_den = np.multiply(inj_den, data_mask)
        inj_frac = np.multiply(inj_frac, data_mask)
        proj_den = np.multiply(proj_den, data_mask)

        # centroid (res = 1 so as to keep in voxel space)
        centroid = mca.calculate_injection_centroid(inj_den,
                                                    inj_frac,
                                                    resolution = 1)

        # compute distances
        #print "Computing Distances"
        distance = compute_distance_from_centroid(centroid, structural_mask)

        # find ratio inside mask
        injection_ratio = sum(inj_den[np.where(functional_mask)]) / sum(inj_den.flatten())
        projection_ratio = sum(proj_den[np.where(functional_mask)]) / sum(proj_den.flatten())

        # add to lists
        centroids += [centroid]
        ratios_in += [injection_ratio]
        distances += [distance]
        projections += [ proj_den[np.where(structural_mask)] ]
        projection_ratios += [projection_ratio]
    return centroids, ratios_in, distances, projections, projection_ratios

def compute_distance_from_centroid(centroid, mask):
    
    # get coordinates for all (returns [ [a,b,..,z],[a,b,...,z],[a,b,...,z] ] format)
    coordinates = np.where( mask )

    # but we want in [ [a,b,c],[a,b,c],...,[a,b,c] ] format
    coordinates = np.stack(coordinates, axis = -1)
    
    # compute distances between centroid and all other points in mask
    # (returns array(array(distances)))
    distances = cdist(np.atleast_2d(centroid),coordinates)

    return distances.ravel()
    
def fit_glm(categorical_var, distances, projections):
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
        groups = np.array(categorical_var)

        dummy = sm.categorical(groups, drop=True)
        x = distances[exp]

        # drop reference category
        X = np.column_stack((x, dummy[:,1:]))
        X = sm.add_constant(X, prepend=False)

        # y Use log projection density
        threshold = 1.6e-4
        projections[exp][projections[exp] < threshold] = 0
        epsilon = 0.01
        y = np.log10( projections[exp] + epsilon )

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
hipp = structure_tree.get_structures_by_acronym(['HPF'])[0]
hipp_mask = mcc.get_structure_mask(hipp['id'])[0]
isohipp_mask = iso_mask + hipp_mask

# grab some experiments
ctx_experiments = pd.DataFrame(mcc.get_experiments(cre=False, 
                                       injection_structure_ids=[iso['id']]))
cre_experiments = pd.DataFrame(mcc.get_experiments(cre=['Emx1-IRES-Cre','Rbp4-Cre_KL100'],
                                      injection_structure_ids = [iso['id']]))
ctx_experiments = pd.concat([ctx_experiments, cre_experiments])
hipp_experiments = pd.DataFrame(mcc.get_experiments(cre=False, 
                                       injection_structure_ids=[hipp['id']]))

# DMN maksks
if platform.system() == 'Windows':
    path = r'2019 DMN'
elif platform.system() == 'Darwin':
    path = r'2019 DMN'
outpath = os.path.join(path, 'data_files')
mask_path = os.path.join(path, 'fMRI_masks', 'all ICA components and masks zscore 1')

for num in np.arange(5):
    print(num)
    mask, _ = nrrd.read(os.path.join(
            mask_path,
            'ica_all_05_icasso_iter_1000_comp_{0}_mask_z_1_allen_masked_sym_thresh_2.nrrd'.format(num)))
    if num == 3:
        in_or_out = shrink_to_mask( mask, isohipp_mask )
        hipp_centroids, hipp_ratios_in, hipp_distances, hipp_projections, hipp_proj_ratios = calculate_centroids_and_distances(
            hipp_experiments, mask, isohipp_mask)
        hipp_d_coeff, hipp_mask_coeff, hipp_tvals, hipp_pvals  = fit_glm(in_or_out, hipp_distances, hipp_projections)
        hipp_glm_dat = pd.DataFrame({'id': hipp_experiments['id'],
                                 'injection structure': hipp_experiments['structure_abbrev'],
                           'injection mask fraction': hipp_ratios_in,
                           'projection mask fraction': hipp_proj_ratios,
                           'distance coefficient': hipp_d_coeff, 
                           'mask coefficient': hipp_mask_coeff,
                           't values': hipp_tvals,
                           'p values': hipp_pvals})
        hipp_glm_dat.to_csv(os.path.join(outpath, 'hippocampal_injections_ICA_{0}.csv'.format(num)))
    # do ctx overlap for all
    ctx_centroids, ctx_ratios_in, ctx_distances, ctx_projections, ctx_proj_ratios = calculate_centroids_and_distances(
        ctx_experiments, mask, iso_mask)
    in_or_out = shrink_to_mask( mask, iso_mask )
    ctx_d_coeff, ctx_mask_coeff, ctx_tvals, ctx_pvals  = fit_glm(in_or_out, ctx_distances, ctx_projections)
    ctx_glm_dat = pd.DataFrame({'id': ctx_experiments['id'],
                                'injection structure': ctx_experiments['structure_abbrev'],
                           'injection mask fraction': ctx_ratios_in,
                           'projection mask fraction': ctx_proj_ratios,
                           'distance coefficient': ctx_d_coeff, 
                           'mask coefficient': ctx_mask_coeff,
                           't values': ctx_tvals,
                           'p values': ctx_pvals})
    ctx_glm_dat.to_csv(os.path.join(outpath, 'ICA_{0}.csv'.format(num)))

        

