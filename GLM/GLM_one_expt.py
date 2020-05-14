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
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi

sns.set_context('paper')
sns.set_style('white')

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

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
        prediction = fit.predict(X)
    return coeff1, coeff2, tvals, pvals, y, prediction
    
# want only inside isocortex
iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
iso_mask = mcc.get_structure_mask(iso['id'])[0]

# grab some experiments
ctx_experiments = pd.DataFrame(mcc.get_experiments(cre=False, 
                                       injection_structure_ids=[iso['id']]))

cre_experiments = pd.DataFrame(mcc.get_experiments(cre=['Emx1-IRES-Cre','Rbp4-Cre_KL100'],
                                      injection_structure_ids = [iso['id']]))
ctx_experiments = pd.concat([ctx_experiments, cre_experiments])

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

# get vector of values inside iso_mask
in_or_out = shrink_to_mask( dmn_mask, iso_mask)
in_or_out_core = shrink_to_mask(core_mask, iso_mask)

# calculate centroids & distances
for exp in [180436360]: #125833030 in #180436360  out

    #print "\n=============    Experiment {}    =============  ".format(exp['id'])
    
    # injection/projection data
    inj_frac = mcc.get_injection_fraction(exp)[0]
    inj_den = mcc.get_injection_density(exp)[0]
    proj_den = mcc.get_projection_density(exp)[0]
    data_mask = mcc.get_data_mask(exp)[0]
    inj_den = np.multiply(inj_den, data_mask)
    
    # centroid (res = 1 so as to keep in voxel space)
    centroid = mca.calculate_injection_centroid(inj_den,
                                                inj_frac,
                                                resolution = 1)
    #print "Centroid is at\t{}".format(centroid)
    
    # compute distances
    #print "Computing Distances"
    distance = compute_distance_from_centroid(centroid,iso_mask)

    proj_iso = proj_den * iso_mask
    
    # find ratio inside mask
    inj_ratio = sum(inj_den[np.where(dmn_mask)]) / sum(inj_den.flatten())
    proj_ratio = sum(proj_iso[np.where(dmn_mask)]) / sum(proj_iso.flatten())
    core_inj_ratio = sum(inj_den[np.where(core_mask)]) / sum(inj_den.flatten())
    core_proj_ratio = sum(proj_iso[np.where(core_mask)]) / sum(proj_iso.flatten())

# fit glm
d_coeff, dmn_coeff, tvals, pvals, yvals, prediction  = fit_glm(in_or_out, [distance], 
                                                        [proj_den[np.where(iso_mask)]] )


#%% Plot
outpath = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/_new_figures/Figure_S2'

fig, ax = plt.subplots(figsize = (1.5, 1.5))
#in_or_out = bool(in_or_out)
sns.scatterplot(distance[np.logical_not(in_or_out)],
               yvals[np.logical_not(in_or_out)], color = 'lime',
               s = 2, alpha = 0.7, ax = ax, rasterized = True, linewidth = 0.5)
sns.scatterplot(distance[np.where(in_or_out)], yvals[np.where(in_or_out)],
                s=2, color='c', alpha = 0.7, ax = ax, rasterized = True, linewidth = 0.5)
sns.lineplot(distance[np.where(in_or_out)], prediction[np.where(in_or_out)],
             color = 'c', ax = ax, label = 'in-DMN voxels')
sns.lineplot(distance[np.logical_not(in_or_out)], 
                prediction[np.logical_not(in_or_out)], color = 'lime',
                ax = ax, label = 'out-DMN voxels')
sns.despine()
ax.tick_params(top = False, right = False, left = True, bottom=True)
ax.set_ylim([-2, 0])
ax.set_xlim([0, 100])
plt.xticks(fontsize = 7)
plt.yticks(fontsize = 7)
ax.set_xlabel('Distance From Centroid (voxels)', fontsize = 7)
ax.set_ylabel('Log10 Projection Density', fontsize = 7)
plt.legend(loc = 1, fontsize = 5)

#sns.scatterplot(distance, in_or_out)
#everything on zorder -1 or lower will be rasterized
plt.savefig(os.path.join(outpath, 'wt_out_GLM_example.pdf'), type = 'pdf',
            bbox_inches='tight', transparent=True, dpi=300)