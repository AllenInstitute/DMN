# -*- coding: utf-8 -*-
"""
Created on Wed Aug 02 16:18:34 2017

@author: jenniferwh
"""
import nrrd
import os
import numpy as np
from scipy.spatial.distance import cdist
import statsmodels.api as sm
from statsmodels.graphics.api import abline_plot
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
import matplotlib.pyplot as plt

mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                            resolution=100)
mca = MouseConnectivityApi()
structure_tree = mcc.get_structure_tree()

def shrink_mask2(mask1, mask2):
    # get coordinates for all (returns [ [a,b,..,z],[a,b,...,z],[a,b,...,z] ] format)
    coordinates = np.where( mask1 )

    # but we want in [ [a,b,c],[a,b,c],...,[a,b,c] ] format
    coordinates = np.stack(coordinates, axis = -1)

    # find mask2 values in mask1 space
    mask_values = []
    for point in coordinates:
        mask_values += [ mask2[ point[0], point[1], point[2] ] ]
    
    return mask_values
    
def compute_distance_from_centroid(centroid, mask):
    
    # get coordinates for all (returns [ [a,b,..,z],[a,b,...,z],[a,b,...,z] ] format)
    coordinates = np.where( mask )

    # but we want in [ [a,b,c],[a,b,c],...,[a,b,c] ] format
    coordinates = np.stack(coordinates, axis = -1)
    
    # compute distances between centroid and all other points in mask
    # (returns array(array(distances)))
    distances = cdist(np.atleast_2d(centroid),coordinates)

    return distances.ravel()
    
    
# want only inside isocortex
iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
iso_mask = mcc.get_structure_mask(iso['id'])[0]

# grab some experiments
ctx_experiments = mcc.get_experiments(cre=False, 
                                       injection_structure_ids=[iso['id']])
dmn_path = r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN paper\fMRI data from Alessandro\fMRI_ICA_map'
dmn_mask, _ = nrrd.read(os.path.join(dmn_path, 'dmn_mask_z_1_allen_masked_sym.nrrd'))

# get vector of values inside iso_mask
in_or_out = shrink_mask2(iso_mask, dmn_mask)

# calculate centroids & distances
centroids = []
ratios_in = []
distances = []
projections = []
for exp in ctx_experiments:

    #print "\n=============    Experiment {}    =============  ".format(exp['id'])
    
    # injection/projection data
    inj_frac = mcc.get_injection_fraction(exp['id'])[0]
    inj_den = mcc.get_injection_density(exp['id'])[0]
    proj_den = mcc.get_projection_density(exp['id'])[0]

    # centroid (res = 1 so as to keep in voxel space)
    centroid = mca.calculate_injection_centroid(inj_den,
                                                inj_frac,
                                                resolution = 1)
    #print "Centroid is at\t{}".format(centroid)
    
    # compute distances
    #print "Computing Distances"
    distance = compute_distance_from_centroid(centroid,iso_mask)
    
    # find ratio inside mask
    ratio = sum(inj_den[np.where(dmn_mask)]) / sum(inj_den.flatten())

    # add to lists
    centroids += [centroid]
    ratios_in += [ratio]
    distances += [distance]
    projections += [ proj_den[np.where(iso_mask)] ]

log_proj = np.log10(projections+np.ones_like(projections))

# fit glm for each experiment
glm_coeff = []
d_coeff = []
for exp in range(len(ctx_experiments)):
    groups = np.array(in_or_out)

    dummy = sm.categorical(groups, drop=True)
    x = distances[exp]

    # drop reference category
    X = np.column_stack((x, dummy[:,1:]))
    #X = sm.add_constant(X, prepend=False)

    # y
    y = log_proj[exp]
    
    # fit
    fit = sm.OLS(y, X).fit()
    
    # add coeff
    glm_coeff += [fit.params[1]]
    d_coeff += [fit.params[0]]

x=ratios_in
y=glm_coeff
dmn_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(fit.summary())

x=ratios_in
y=d_coeff
distance_fit = sm.OLS(y, sm.add_constant(x, prepend=True)).fit()
print(fit.summary())

fig, ax = plt.subplots(2,1, sharex = True, figsize = (10,10))
fig.suptitle('DMN projections', y = 0.93, fontsize = 18)
ax[0].scatter(ratios_in, glm_coeff, c = 'r')
abline_plot(model_results=dmn_fit, ax = ax[0], c = 'k')
ax[1].set_ylabel('Distance Coefficient (Log distance)')
ax[0].set_ylabel('DMN Coefficient')

ax[1].scatter(ratios_in, d_coeff)
abline_plot(model_results=distance_fit, ax = ax[1], c = 'g')
plt.xlabel('Fraction of injection in DMN')
plt.show()