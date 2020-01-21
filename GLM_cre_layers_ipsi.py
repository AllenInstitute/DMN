#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:39:34 2020

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

if platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN'
elif platform.system() == 'Darwin':
    path = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN'
    
maskpath = os.path.join(path, 'fMRI_masks')
datpath = os.path.join(path, 'data_files')

dat = pd.read_csv(os.path.join(datpath, 'Ext Data Table 3_Corticocortical NPV matrices.csv'))
dat.rename(columns = {'Exp image series ID': 'image_series_id', 'Rbp4 anchor image series ID': 'Rbp4_anchor_id',
                     'distance from anchor': 'distance', 'Consensus Rbp4 anchor Source': 'Rbp4_source'}, 
           inplace = True)
print(len(dat))
print(len(dat['Rbp4_anchor_id'].unique()))

dat.loc[dat['Mouse Line'] == 'C57BL/6J / Emx1', 'layer'] = 'all'
dat.loc[dat['Mouse Line'].isin(['Cux2-IRES-Cre', 'Sepw1-Cre_NP39']), 'layer'] = '2/3'
dat.loc[dat['Mouse Line'].isin(['Nr5a1-Cre', 'Scnn1a-Tg3-Cre', 'Rorb-IRES2-Cre']), 'layer'] = '4'
dat.loc[dat['Mouse Line'].isin(['Tlx3-Cre_PL56']), 'layer'] = '5 IT' #separate
dat.loc[dat['Mouse Line'].isin(['Rbp4-Cre_KL100']), 'layer'] = '5 IT PT' #separate
dat.loc[dat['Mouse Line'].isin(['Chrna2-Cre_OE25', 'Efr3a-Cre_NO108', 'Sim1-Cre_KJ18',
                               'A930038C07Rik-Tg1-Cre']), 'layer'] = '5 PT'
dat.loc[dat['Mouse Line'].isin(['Ntsr1-Cre_GN220', 'Syt6-Cre_KI148']), 'layer'] = '6'

masks, _ = nrrd.read(os.path.join(maskpath, 'dmn_mask_and_core.nrrd'))
dmn_mask = np.zeros(masks.shape)
dmn_mask[np.where(masks > 0)] = 1
dmn_mask = dmn_mask[:,:,:57]
core_mask = np.zeros(masks.shape)
core_mask[np.where(masks == 2)] = 1
core_mask = core_mask[:,:,:57]

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
        y = np.log10 ( projections[exp] + 3.1e-14 ) #1/2 min proj. value
    
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
iso_mask = iso_mask[:,:,:57]

ctx_experiments = pd.DataFrame(mcc.get_experiments(cre=False, 
                                       injection_structure_ids=[iso['id']]))
mouselines = dat['Mouse Line'].unique()
mouselines = np.append(mouselines, 'Emx1-IRES-Cre')
cre_experiments = pd.DataFrame(mcc.get_experiments(cre=mouselines,
                                      injection_structure_ids = [iso['id']]))
cre_experiments = cre_experiments[cre_experiments['id'].isin(
    dat['image_series_id'].unique())]
ctx_experiments = pd.concat([ctx_experiments, cre_experiments])
print(len(ctx_experiments))

# get vector of values inside iso_mask
in_or_out = shrink_mask2(iso_mask, dmn_mask)
in_or_out_core = shrink_mask2(iso_mask, core_mask)

# calculate centroids & distances
centroids = []
inj_ratios = []
proj_ratios = []
core_inj_ratios = []
core_proj_ratios = []
distances = []
projections = []
for exp in ctx_experiments['id']:

    #print "\n=============    Experiment {}    =============  ".format(exp['id'])
    
    # injection/projection data
    inj_frac = mcc.get_injection_fraction(exp)[0]
    inj_den = mcc.get_injection_density(exp)[0]
    proj_den = mcc.get_projection_density(exp)[0]
    
    # centroid (res = 1 so as to keep in voxel space)
    centroid = mca.calculate_injection_centroid(inj_den,
                                                inj_frac,
                                                resolution = 1)
    if centroid[2] > 57: #right side injection
        inj_frac = np.flip(inj_frac, 2)
        inj_den = np.flip(inj_den, 2)
        proj_den = np.flip(proj_den, 2)
    inj_frac = inj_frac[:,:,:57]
    inj_den = inj_den[:,:,:57]
    proj_den = proj_den[:,:,:57]
    
    # compute distances
    #print "Computing Distances"
    distance = compute_distance_from_centroid(centroid,iso_mask)
    
    inj_iso = inj_den * iso_mask
    proj_iso = proj_den * iso_mask
    
    # find ratio inside mask
    inj_ratio = sum(inj_iso[np.where(dmn_mask)]) / sum(inj_iso.flatten())
    proj_ratio = sum(proj_iso[np.where(dmn_mask)]) / sum(proj_iso.flatten())
    core_inj_ratio = sum(inj_iso[np.where(core_mask)]) / sum(inj_iso.flatten())
    core_proj_ratio = sum(proj_iso[np.where(core_mask)]) / sum(proj_iso.flatten())

    # add to lists
    centroids += [centroid]
    inj_ratios += [inj_ratio]
    proj_ratios += [proj_ratio]
    core_inj_ratios += [core_inj_ratio]
    core_proj_ratios += [core_proj_ratio]
    distances += [distance]
    projections += [ proj_den[np.where(iso_mask)] ]

# fit glm
d_coeff, dmn_coeff, tvals, pvals  = fit_glm(in_or_out, distances, projections)
d_coeff_core, dmn_coeff_core, tvals_core, pvals_core  = fit_glm(in_or_out_core, distances, projections)

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

ctx_glm_dat = pd.DataFrame({'id': ctx_experiments['id'], 
                       'injection structure': ctx_experiments['structure_abbrev'],
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
ctx_glm_dat.to_csv(os.path.join(datpath, 'wt_cre_matched_layer_injections_DMN_and_core_projections_coefficients_ipsi.csv'),
              index = False)