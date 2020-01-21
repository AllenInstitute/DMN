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

# grab some experiments
ctx_experiments = pd.DataFrame(mcc.get_experiments(cre=False, 
                                       injection_structure_ids=[iso['id']]))

cre_experiments = pd.DataFrame(mcc.get_experiments(cre=['Emx1-IRES-Cre','Rbp4-Cre_KL100'],
                                      injection_structure_ids = [iso['id']]))
ctx_experiments = pd.concat([ctx_experiments, cre_experiments])

fail_expts = [114008926, 120280939, 180073473, 180403712, 180601025, 183174303, 183329222,
              249396394, 296047806, 299446445, 301060890, 303784745, 480069939, 482578964, 
              506947040, 514333422, 525796603, 545428296, 559878074, 638314843, 182888003,
             304585910, 183171679, 272930013, 523718075, 517072832, 148964212, 304762965,
             566992832, 272930013, 304762965, 266250904, 114399224, 286483411, 286417464,
             593277684, 546103149, 642809043, 286483411, 304564721] #VISp outlier excluded

ctx_experiments = ctx_experiments[~ctx_experiments['id'].isin(fail_expts)]

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
    #print "Centroid is at\t{}".format(centroid)
    
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
ctx_glm_dat.to_csv(os.path.join(outpath, 'wt_cre_ctx_injections_DMN_and_core_projections_coefficients.csv'),
              index = False)

#%% PLot

pvals = [dmn_fit.pvalues[1],
         core_fit.pvalues[1],
         distance_fit.pvalues[1]]
fdrcorr = sm.stats.fdrcorrection(pvals, alpha=0.05, method='indep')
print(fdrcorr)

ctx_glm_dat['inj_percent_dmn'] = ctx_glm_dat['injection dmn fraction'] * 100
ctx_glm_dat['proj_percent_dmn'] = ctx_glm_dat['projection dmn fraction'] * 100
ctx_glm_dat['inj_percent_core'] = ctx_glm_dat['injection core fraction'] * 100
ctx_glm_dat['proj_percent_core'] = ctx_glm_dat['projection core fraction'] * 100
fig, ax = plt.subplots(figsize = (3, 3))
sns.regplot(x = 'inj_percent_dmn',
            y = 'DMN coefficient',
            data = ctx_glm_dat,
            ax = ax,
            color = 'gray',
            scatter_kws={'alpha':0.7, 's':10},
           label = r'DMN ($R^2$={0})'.format(round(dmn_fit.rsquared, 2)))
sns.regplot(x = 'inj_percent_core',
            y = 'DMN core coefficient',
            data = ctx_glm_dat,
            ax = ax,
            color = 'k',
            scatter_kws={'alpha':0.7, 's':10},
           label = r'Core ($R^2$={0})'.format(round(core_fit.rsquared, 2)))
sns.regplot(x = 'inj_percent_core',
            y = 'distance coefficient',
            data = ctx_glm_dat,
            ax = ax,
            color = 'w',
            scatter_kws={'edgecolor': 'k', 's':10},
            line_kws={'color':'k', 'zorder': -1, 'alpha':0.8},
           label = r'Distance ($R^2$={0})'.format(round(distance_fit.rsquared, 2)))
ax.set_ylabel('Coefficient', fontsize = 8)
ax.yaxis.labelpad = 1
ax.xaxis.labelpad = 1
ax.set_xlabel('Injection % DMN', fontsize = 8)

leg = plt.legend(fontsize = 6, labelspacing=-0.1, frameon = False,
          bbox_to_anchor = [0.5, 0.8])
for lh in leg.legendHandles: 
    lh.set_alpha(1)
plt.xticks([0, 50, 100], fontsize = 8)
sns.despine()
ax.tick_params(top = False, right = False, pad = -2)

#everything on zorder -1 or lower will be rasterized
ax.set_rasterization_zorder(0)

#plt.savefig(os.path.join(path, 'wt_distance_DMN_coeff.pdf'), bbox_inches='tight', transparent=True, dpi=1000)