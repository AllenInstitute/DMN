#!/usr/bin/env python
# coding: utf-8

# In[1]:
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns
sns.set_context('paper')
sns.set_style('white')

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
# In[2]:
import platform
if platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\Figure 4 DMN coefficient wt'
    maskpath = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\fMRI masks'
elif platform.system() == 'Darwin':
    path = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/Figure 4 DMN coefficient wt'
    maskpath = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/fMRI masks'
# In[3]:
import nrrd
import os
both_dmn_masks, _ = nrrd.read(os.path.join(maskpath, 'dmn_mask_and_core.nrrd'))
dmn_mask = np.zeros(both_dmn_masks.shape)
dmn_mask[np.where(both_dmn_masks > 0)] = 1
core_mask = np.zeros(both_dmn_masks.shape)
core_mask[np.where(both_dmn_masks == 2)] = 1

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json',
                            resolution=100)
# 100 um for calculations
structure_tree = mcc.get_structure_tree()
iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0] # get a 100 um isocortex mask
iso_mask = mcc.get_structure_mask(iso['id'])[0]


# In[19]:

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


# In[20]:

# wild type experiments
ctx_experiments = mcc.get_experiments(cre=False, 
                                       injection_structure_ids=[iso['id']])
print(len(ctx_experiments))

# In[21]:
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
mca = MouseConnectivityApi()

def calculate_centroids_and_distances(experiments_list, structure_mask, functional_mask):
    centroids = []
    injection_ratios_in = []
    projection_ratios_in = []
    distances = []
    projections = []
    for exp in experiments_list:
        # injection/projection data
        inj_frac = mcc.get_injection_fraction(exp['id'])[0]
        inj_den = mcc.get_injection_density(exp['id'])[0]
        proj_den = mcc.get_projection_density(exp['id'])[0]
        data_mask = mcc.get_data_mask(exp['id'])[0]

        inj_den = np.multiply(inj_den, data_mask)
        inj_frac = np.multiply(inj_frac, data_mask)
        proj_den = np.multiply(proj_den, data_mask)

        # compute normalized projection density
        npd = proj_den / inj_den.sum()

        # centroid (res = 1 so as to keep in voxel space)
        centroid = mca.calculate_injection_centroid(inj_den,
                                                    inj_frac,
                                                    resolution = 1)

        # compute distances
        distance = compute_distances_from_centroid(centroid, structure_mask)

        # find ratio inside mask
        inj_ratio = sum(inj_den[np.where(functional_mask)]) / sum(inj_den.flatten())
        ctx_proj = np.multiply(proj_den, structure_mask)
        proj_ratio = sum(ctx_proj[np.where(functional_mask)]) / sum(ctx_proj.flatten())

        # add to lists
        centroids += [centroid]
        injection_ratios_in += [inj_ratio]
        projection_ratios_in += [proj_ratio]
        distances += [distance]
        projections += [ npd[np.where(structure_mask)] ]
    
    return centroids, injection_ratios_in, projection_ratios_in, distances, projections


# In[22]:
ctx_centroids, ctx_injection_ratios_in, ctx_projection_ratios_in, ctx_distances, ctx_projections = calculate_centroids_and_distances(ctx_experiments, iso_mask, dmn_mask)
ctx_core_centroids, ctx_core_injection_ratios_in, ctx_core_projection_ratios_in, ctx_core_distances, ctx_core_projections = calculate_centroids_and_distances(ctx_experiments, iso_mask, core_mask)

# In[23]:
import statsmodels.api as sm
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
        y = np.log10 ( projections[exp] + 0.0000001 )
    
        # fit
        fit = sm.OLS(y, X).fit()
    
        # add coeff
        coeff1 += [fit.params[0]]
        coeff2 += [fit.params[1]]
        tvals += [fit.tvalues]
        pvals += [fit.pvalues]
    return coeff1, coeff2, tvals, pvals


# In[24]:
# get vector of dmn values inside iso_mask
in_or_out_ctx = shrink_to_mask(dmn_mask, iso_mask)
in_or_out_core_ctx = shrink_to_mask(core_mask, iso_mask)


# In[25]:
ctx_d_coeff, ctx_dmn_coeff, ctx_tvals, ctx_pvals  = fit_glm(in_or_out_ctx, ctx_distances, ctx_projections)
ctx_core_d_coeff, ctx_core_coeff, ctx_core_tvals, ctx_core_pvals = fit_glm(in_or_out_core_ctx, ctx_core_distances, ctx_core_projections)

# In[26]:
import pandas as pd
ctx_experiments = pd.DataFrame(ctx_experiments)

# In[28]:
ctx_glm_dat = pd.DataFrame({'id': ctx_experiments['id'],
                           'injection structure': ctx_experiments['structure_abbrev'],
                           'injection dmn fraction': ctx_injection_ratios_in,
                           'projection dmn fraction': ctx_projection_ratios_in,
                           'distance coefficient': ctx_d_coeff, 
                           'DMN coefficient': ctx_dmn_coeff,
                           'DMN t values': ctx_tvals,
                           'DMN p values': ctx_pvals,
                           'injection core fraction': ctx_core_injection_ratios_in,
                           'projection core fraction': ctx_core_projection_ratios_in,
                           'core distance coefficient': ctx_core_d_coeff,
                           'DMN core coefficient': ctx_core_coeff,
                           'core t values': ctx_core_tvals,
                           'core p values': ctx_core_pvals})


# In[29]:
ctx_glm_dat.to_csv(
        os.path.join(path, 
                     'wt_ctx_injections_DMN_and_core_projections_coefficients.csv'),
                     index = False)


# In[30]:
ctx_dmn_fit = sm.OLS(ctx_dmn_coeff, sm.add_constant(ctx_injection_ratios_in, prepend=True)).fit()
ctx_d_fit = sm.OLS(ctx_d_coeff, sm.add_constant(ctx_injection_ratios_in, prepend=True)).fit()
# In[31]:
print(ctx_dmn_fit.summary())
# In[32]:
print(ctx_d_fit.summary())

# In[33]
core_fit = sm.OLS(ctx_core_coeff, sm.add_constant(ctx_core_injection_ratios_in, prepend=True)).fit()
core_d_fit = sm.OLS(ctx_core_d_coeff, sm.add_constant(ctx_core_injection_ratios_in, prepend=True)).fit()
# In[31]:
print(core_fit.summary())
# In[32]:
print(core_d_fit.summary())
# In[33]:

from allensdk.api.queries.ontologies_api import OntologiesApi
oapi = OntologiesApi()
summary_structures = oapi.get_structures(structure_set_names="'Brain – Summary Structures'")
summary_structure_ids = [item['id'] for item in summary_structures]
print(len(summary_structure_ids))
ctx_exps = pd.DataFrame(mcc.get_experiments(cre=False, injection_structure_ids = [iso['id']]))
print(len(ctx_exps))


# In[36]:


ctx_exps.tail()


# In[37]:


ia_map = structure_tree.get_id_acronym_map()
ai_map = {value: key for key, value in ia_map.iteritems()}


# In[38]:


inj_structure_abbrevs = []
for experiment in ctx_exps['injection_structures'].values:
    abbrevs = [ai_map[structure] for structure in experiment]
    inj_str_abbrev = ''
    for abbrev in abbrevs:
        inj_str_abbrev = inj_str_abbrev + '|' + abbrev
    inj_structure_abbrevs.append(inj_str_abbrev)
ctx_exps['injection structures abbreviations'] = inj_structure_abbrevs


# In[39]:


# Get the value for dmn fraction and DMN coefficient from the dataframe calculated at 100 um for each experiment
percent_dmn = []
dmn_coeff = []
for isid in ctx_exps['id']:
    percent_dmn.append(ctx_glm_dat[ctx_glm_dat['id'] == isid]['projection dmn fraction'].values[0] * 100)
    dmn_coeff.append(ctx_glm_dat[ctx_glm_dat['id'] == isid]['DMN coefficient'].values[0])
ctx_exps['percent_dmn'] = percent_dmn
ctx_exps['DMN_coefficient'] = dmn_coeff


# In[40]:


ctx_exps.keys()


# In[42]:


ctx_exps.head()


# In[43]:


def get_injection_volume(df):
    for experiment in df['id']:
        unionize = mcc.get_experiment_structure_unionizes(experiment_id = experiment, 
                                           hemisphere_ids = [3],
                                           is_injection = True)
        inj_size = unionize[unionize['structure_id'] == 997]['projection_volume'].values[0]
        df.loc[df['id'] == experiment, 'injection_volume'] = inj_size
    return df


# In[44]:


ctx_exps = get_injection_volume(ctx_exps)


# In[45]:


template_ax = np.percentile(template[:,:120,:], axis = 1, q=100)
template_sag = np.percentile(template[:,:,:170], axis = 2, q=100)


# ## Figure 3b

# In[46]:


fig,ax = plt.subplots(1,1)
mask = get_mask('horizontal')
ax.imshow(template_ax, cmap='gray')
ax.imshow(mask, cmap = 'gist_gray', vmax = 3, alpha = 0.8)
ax.imshow(template_ax, cmap = 'gray', alpha = 0.6)
cax = ax.scatter(ctx_exps.injection_z/25,
           ctx_exps.injection_x/25,
           c = ctx_exps['percent_dmn'],
           edgecolor = 'none',
           s = ctx_exps.injection_volume.values*500,
           cmap = 'gray',
                 alpha = 0.9)

cbar = fig.colorbar(cax, ax = ax, orientation='vertical', 
                    fraction=0.046, pad=0.04)
cbar_ax = fig.axes[-1]
cbar.solids.set_rasterized(True)
cbar.solids.set_edgecolor("face")
cbar.set_alpha(1)
cbar.draw_all()
ax.set_aspect('equal')
plt.axis('off')
fig.set_size_inches(10,10)
plt.savefig(os.path.join(path, 'wt_injections_on_DMN_horizontal_cmap_gray.pdf'), 
            bbox_inches='tight', pad_inches=0.3, format='pdf', transparent=True, dpi=1000)


# In[58]:


fig,ax = plt.subplots(1,1)
mask = get_mask('horizontal')
ax.imshow(template_ax, cmap='gray')
ax.imshow(mask, cmap = 'gist_gray', vmax = 3, alpha = 0.8)
ax.imshow(template_ax, cmap = 'gray', alpha = 0.6)
cax = ax.scatter(ctx_exps.injection_z/25,
           ctx_exps.injection_x/25,
           c = ctx_exps['DMN_coefficient'],
           edgecolor = 'none',
           s = ctx_exps.injection_volume.values*50,
           cmap = 'bwr',
                 vmax = 2,
                 vmin = -2,
                 alpha = 0.9)

cbar = fig.colorbar(cax, ax = ax, orientation='vertical', 
                    fraction=0.046, pad=0.04)
cbar_ax = fig.axes[-1]
cbar.solids.set_rasterized(True)
cbar.solids.set_edgecolor("face")
cbar.set_alpha(1)
cbar.draw_all()
ax.set_aspect('equal')
plt.axis('off')
fig.set_size_inches(4,4)
plt.savefig(os.path.join(path, 'wt_injections_on_DMN_horizontal_DMN_coefficient_cmap_bwr.pdf'), 
            bbox_inches='tight', pad_inches=0.3, format='pdf', transparent=True, dpi=300)


# In[48]:


fig,ax = plt.subplots(1,1)
mask = get_mask('sagittal')
ax.imshow(template_sag, cmap = 'gray')
ax.imshow(mask, cmap = 'gist_gray', vmax = 3, alpha = 0.8)
ax.imshow(template_sag, cmap = 'gray', alpha = 0.6)
cax = ax.scatter(ctx_exps.injection_y/25,
           ctx_exps.injection_x/25,
           c = ctx_exps['percent_dmn'],
            cmap = 'gray',
           edgecolor = 'none',
           s = ctx_exps.injection_volume.values*500,
          alpha = 0.8)
cbar = fig.colorbar(cax, ax = ax, orientation='vertical', 
                    fraction=0.046, pad=0.04)
cbar_ax = fig.axes[-1]
cbar.solids.set_rasterized(True)
cbar.solids.set_edgecolor("face")
cbar.set_alpha(1)
cbar.draw_all()
ax.set_aspect('equal')
plt.axis('off')
fig.set_size_inches(10, 10)
plt.savefig(os.path.join(path, 'wt_dmn_proj_injection_spatial_map_sag_cmap_gray.pdf'), 
           bbox_inches='tight', pad_inches=0.3, format='pdf', transparent = True, dpi=1000)


# In[59]:


fig,ax = plt.subplots(1,1)
mask = get_mask('sagittal')
ax.imshow(template_sag, cmap = 'gray')
ax.imshow(mask, cmap = 'gist_gray', vmax = 3, alpha = 0.8)
ax.imshow(template_sag, cmap = 'gray', alpha = 0.6)
cax = ax.scatter(ctx_exps.injection_y/25,
           ctx_exps.injection_x/25,
           c = ctx_exps['DMN_coefficient'],
            cmap = 'bwr',
           edgecolor = 'none',
           s = ctx_exps.injection_volume.values*50,
          alpha = 0.8)
cbar = fig.colorbar(cax, ax = ax, orientation='vertical', 
                    fraction=0.046, pad=0.04)
cbar_ax = fig.axes[-1]
cbar.solids.set_rasterized(True)
cbar.solids.set_edgecolor("face")
cbar.set_alpha(1)
cbar.draw_all()
ax.set_aspect('equal')
plt.axis('off')
fig.set_size_inches(4,4)
plt.savefig(os.path.join(path, 'wt_dmn_proj_injection_spatial_map_sag_DMN_coefficient_cmap_RdGy_r.pdf'), 
           bbox_inches='tight', pad_inches=0.3, format='pdf', transparent = True, dpi=300)


# ## Boxplots

# ## Isocortex Injections

# In[50]:


ctx_glm_dat = ctx_glm_dat.merge(ctx_exps[['id', 'injection structures abbreviations']], on='id')


# In[51]:


ctx_glm_dat.head()


# In[52]:


strs = mcc.rank_structures(ctx_glm_dat['id'], is_injection = True)
strids = [exp[0]['structure_id'] for exp in strs]


# In[168]:


sorted_names = pd.DataFrame(
    ctx_glm_dat[['DMN coefficient', 
        'injection structure', 
        'injection dmn fraction', 
        'projection dmn fraction']].groupby(
        'injection structure').mean().sort_values(
        by = 'DMN coefficient', ascending = False))


# In[169]:


sorted_names.index


# In[170]:


sd = ctx_glm_dat[['DMN coefficient', 
        'injection structure']].groupby(
        'injection structure').std().reset_index()
sd.rename(columns = {'DMN coefficient': 'DMN coefficient STD'}, inplace = True)
sorted_names = sorted_names.merge(sd, on='injection structure', how = 'left')


# In[171]:


count = ctx_glm_dat[['DMN coefficient', 
        'injection structure']].groupby(
        'injection structure').count().reset_index()
count.rename(columns = {'DMN coefficient': 'Number of Experiments'}, inplace = True)
sorted_names = sorted_names.merge(count, on='injection structure', how = 'left')


# In[172]:


sorted_names.head()


# In[173]:


sorted_names.to_csv(os.path.join(path, 'wt mean dmn metrics by structure.csv'), index = False)


# In[174]:


sorted_names = pd.DataFrame(
    ctx_glm_dat[['DMN coefficient', 
        'injection structure', 
        'injection dmn fraction', 
        'projection dmn fraction']].groupby(
        'injection structure').mean().sort_values(
        by = 'DMN coefficient', ascending = False))
sorted_names.to_csv(os.path.join(path, 'wt experiment count by structure.csv'))


# In[175]:


mcd = {'red': '#ff0000', 'orange': '#f9922b', 'yellow': '#ffff66', 'light_blue': '#90bff9',
      'dark_blue': '#5252a9', 'purple': '#7c429b'}


# In[176]:


mod_colors = [mcd['red'], mcd['red'], mcd['red'], mcd['red'], mcd['red'], mcd['dark_blue'], mcd['red'],
             mcd['dark_blue'], mcd['orange'], mcd['dark_blue'], mcd['light_blue'], mcd['orange'],
             mcd['light_blue'], mcd['light_blue'], mcd['dark_blue'], mcd['orange'], mcd['purple'],
             mcd['light_blue'], mcd['orange'], mcd['light_blue'], mcd['dark_blue'], mcd['yellow'],
             mcd['purple'], mcd['orange'], mcd['red'], mcd['purple'], mcd['light_blue'], mcd['purple'],
             mcd['yellow'], mcd['orange'], mcd['yellow'], mcd['yellow'], mcd['orange'], mcd['orange'],
             mcd['orange'], mcd['red'], mcd['orange']]


# In[177]:


sorted_names.index


# In[178]:


sns.set(font_scale=2)
sns.set_style('white')


# In[179]:


fig, ax = plt.subplots(1, figsize = (15, 5))
sns.boxplot('injection structure', 'DMN coefficient', data = ctx_glm_dat, 
            order = sorted_names.index, ax = ax, color = 'white')
sns.stripplot('injection structure', 'DMN coefficient', data = ctx_glm_dat, 
            order = sorted_names.index, ax = ax, palette = mod_colors, s = 8, alpha = 0.7)
ax.set_xticks(np.linspace(ax.get_xbound()[0]+0.5, ax.get_xbound()[1]-0.5, len(sorted_names.index)))
plt.xticks(rotation = -90)
sns.despine()
ax.tick_params(top = False, right = False)
ax.set_ylabel("")
ax.set_xlabel("")
plt.title('Cortical Structure', y=1.03)
ax.axhline(y=0, xmin=0, xmax=100, color = 'grey', linestyle = 'dashed', zorder = -1)
plt.savefig(os.path.join(path, 'DMN_proj_by_source_wt_boxplot_horizontal.pdf'), 
            bbox_inches='tight', pad_inches=0.3, format='pdf', transparent=True, dpi=1000)


# In[144]:


ctx_glm_dat.loc[ctx_glm_dat['injection structure'].isin(['ACAd', 'ACAv', 'ORBvl', 'ORBm', 'ILA', 
                                                    'PL', 'ORBl', 'FRP']), 'module'] = 'Prefrontal'
ctx_glm_dat.loc[ctx_glm_dat['injection structure'].isin(['VISam', 'RSPagl', 'RSPv', 'VISpm',
                                                        'RSPd']), 'module'] = 'Medial'
ctx_glm_dat.loc[ctx_glm_dat['injection structure'].isin(['SSp-tr', 'SSp-ll', 'SSp-bfd', 'SSp-un',
                                                        'SSp-ul', 'MOs', 'MOp', 'SSs', 'SSp-n',
                                                        'SSp-m']), 'module'] = 'Somatomotor'
ctx_glm_dat.loc[ctx_glm_dat['injection structure'].isin(['VISal', 'VISl', 'VISli', 'VISp', 'VISpl',
                                                        'VISpor']), 'module'] = 'Visual'
ctx_glm_dat.loc[ctx_glm_dat['injection structure'].isin(['AUDpo', 'AUDd', 'AUDp', 'AUDv']), 
                'module'] = 'Auditory'
ctx_glm_dat.loc[ctx_glm_dat['injection structure'].isin(['ECT', 'GU', 'AId', 'VISC']), 
                'module'] = 'Lateral'


# In[145]:


mod_order = ['Prefrontal', 'Medial', 'Visual', 'Auditory', 'Somatomotor', 'Lateral', ]
colors = [mcd['red'], mcd['dark_blue'], mcd['light_blue'], mcd['purple'], mcd['orange'], mcd['yellow']]


# In[180]:


fig, ax = plt.subplots(1, figsize = (4, 5))
sns.boxplot('module', 'DMN coefficient', data = ctx_glm_dat, 
            order = mod_order, ax = ax, color = 'white')
sns.stripplot('module', 'DMN coefficient', data = ctx_glm_dat, 
            order = mod_order, ax = ax, palette = colors, s = 8, alpha = 0.7)
sns.despine()
ax.tick_params(top = False, right = False)
plt.xticks(rotation = -90)
plt.title('Cortical Module', y=1.03)
ax.set_xlabel("")
ax.axhline(y=0, xmin=0, xmax=100, color = 'grey', linestyle = 'dashed', zorder = -1)
plt.savefig(os.path.join(path, 'DMN_proj_by_module_wt_boxplot_horizontal.pdf'), 
            bbox_inches='tight', pad_inches=0.3, format='pdf', transparent=True, dpi=1000)


# In[114]:


fig, ax = plt.subplots(1, figsize = (2, 7))
sns.boxplot('DMN coefficient', 'injection structure', data = ctx_glm_dat, order = sorted_names.index, 
            ax = ax, color = 'white')
sns.stripplot('DMN coefficient', 'injection structure', data = ctx_glm_dat, order = sorted_names.index, 
            ax = ax, palette = mod_colors, s = 5, alpha = 0.7)
ax.set_yticks(np.linspace(ax.get_ybound()[0]+0.5, ax.get_ybound()[1]-0.5, len(sorted_names.index)))
plt.title('DMN Projection Strength by Injected Source')
ax.axvline(x=0, ymin=0, ymax=100, color = 'grey', linestyle = 'dashed', zorder = -1)
plt.savefig(os.path.join(path, 'DMN_proj_by_source_wt_boxplot.pdf'), 
           bbox_inches='tight', pad_inches=0.3, format='pdf', transparent = True, dpi=1000)


# In[115]:


sorted_names = pd.DataFrame(
    ctx_glm_dat[['DMN coefficient', 
        'injection structures abbreviations', 
        'injection dmn fraction', 
        'projection dmn fraction']].groupby(
        'injection structures abbreviations').mean().sort_values(
        by = 'DMN coefficient', ascending = False))


# In[116]:


fig, ax = plt.subplots(1, figsize = (5, 23))
sns.boxplot('DMN coefficient', 'injection structures abbreviations', data = ctx_glm_dat, order = sorted_names.index, 
            ax = ax, palette = 'inferno_r')
ax.set_yticks(np.linspace(ax.get_ybound()[0]+0.5, ax.get_ybound()[1]-0.5, len(sorted_names.index)))
plt.title('DMN Projection Strength by Injected Source')
ax.axvline(x=0, ymin=0, ymax=100, color = 'grey', linestyle = 'dashed', zorder = -1)
plt.savefig(os.path.join(path, 'DMN_proj_by_multi_source_wt_boxplot.pdf'), 
           bbox_inches='tight', pad_inches=0.3, format='pdf', transparent = True, dpi=1000)


# ## Thalamus injections

# In[ ]:


mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                            resolution=100) # Switch to 100 um for calculations


# In[ ]:


iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0] # get a 100 um isocortex mask
iso_mask = mcc.get_structure_mask(iso['id'])[0]


# In[ ]:


thal = structure_tree.get_structures_by_acronym(['TH'])[0]
th_mask = mcc.get_structure_mask(thal['id'])[0]
isoth_mask = th_mask + iso_mask
thal_exps = mcc.get_experiments(cre=False, injection_structure_ids = [thal['id']])
print(len(thal_exps))


# In[ ]:


ptlp = structure_tree.get_structures_by_acronym(['PTLp'])[0]
iso_descendants = structure_tree.descendants([iso['id']])[0]
ptlp_descendants = structure_tree.descendants([ptlp['id']])[0]
iso_descendant_ids = [region['id'] for region in iso_descendants]
ptlp_descendant_ids = [region['id'] for region in ptlp_descendants]
descendants = iso_descendant_ids + ptlp_descendant_ids


# In[ ]:


## Eliminate experiments with significant cortical leakage
leakage_experiments = []
for experiment in thal_exps:
    ij_strs = [structure['id'] for structure in experiment['injection-structures']]
    ctx_strs_in_th_injection = [i for i in ij_strs if i in descendants]
    if len(ctx_strs_in_th_injection) > 0:
        leakage_experiments.append(experiment['id'])
print(len(leakage_experiments))


# In[ ]:


thal_exps = [exp for exp in thal_exps if exp['id'] not in leakage_experiments]
print(len(thal_exps))


# In[ ]:


(thal_centroids, thal_injection_ratios_in, 
 thal_projection_ratios_in, thal_distances, 
 thal_projections) = calculate_centroids_and_distances(thal_exps, isoth_mask)


# In[ ]:


in_or_out_th = shrink_to_mask(dmn_mask, isoth_mask)
print(len(isoth_mask))


# In[ ]:


thal_d_coeff, thal_dmn_coeff, thal_tvals, thal_pvals  = fit_glm(in_or_out_th, thal_distances, thal_projections)


# In[ ]:


thal_exps = pd.DataFrame(thal_exps)
th_glm_dat = pd.DataFrame({'id': thal_exps['id'],
                           'injection structure': thal_exps['structure-abbrev'],
                           'injection dmn fraction': thal_injection_ratios_in,
                           'projection dmn fraction': thal_projection_ratios_in,
                           'distance coefficient': thal_d_coeff, 
                           'DMN coefficient': thal_dmn_coeff,
                           't values': thal_tvals,
                           'p values': thal_pvals})


# In[ ]:


percent_dmn = []
dmn_coeff = []
for isid in thal_exps['id']:
    percent_dmn.append(th_glm_dat[th_glm_dat['id'] == isid]['projection dmn fraction'].values[0] * 100)
    dmn_coeff.append(th_glm_dat[th_glm_dat['id'] == isid]['DMN coefficient'].values[0])
thal_exps['percent_dmn'] = percent_dmn
thal_exps['DMN_coefficient'] = dmn_coeff
thal_exps = separate_coordinates(thal_exps)
thal_exps = get_injection_volume(thal_exps)


# In[ ]:


th_glm_dat.to_csv(os.path.join(path, 'thalamic injections DMN distance coefficient ctx and thal DMN mask.csv'))


# In[ ]:


strs = mcc.rank_structures(th_glm_dat['id'], is_injection = True)
strids = [exp[0]['structure_id'] for exp in strs]
sorted_names = pd.DataFrame(
    th_glm_dat[['DMN coefficient', 
        'injection structure', 
        'injection dmn fraction', 
        'projection dmn fraction']].groupby(
        'injection structure').mean().sort_values(
        by = 'DMN coefficient', ascending = False))


# In[ ]:


fig, ax = plt.subplots(1, figsize = (9, 5))
sns.boxplot('injection structure', 'DMN coefficient', data = th_glm_dat, 
            order = sorted_names.index, ax = ax, palette = 'inferno_r')
ax.set_xticks(np.linspace(ax.get_xbound()[0]+0.5, ax.get_xbound()[1]-0.5, len(sorted_names.index)))
plt.xticks(rotation = -90)
plt.title('Thalamus Injections', y=1.03)
ax.axhline(y=0, xmin=0, xmax=100, color = 'grey', linestyle = 'dashed', zorder = -1)
ax.set_ylim(-0.5, 1)
plt.savefig(os.path.join(path, 'DMN_proj_by_source_wt_thalamus_boxplot_horizontal.pdf'), 
            bbox_inches='tight', pad_inches=0.3, format='pdf', transparent=True, dpi=1000)


# In[ ]:


thal_dmn_fit = sm.OLS(thal_dmn_coeff, sm.add_constant(thal_injection_ratios_in, prepend=True)).fit()
thal_d_fit = sm.OLS(thal_d_coeff, sm.add_constant(thal_injection_ratios_in, prepend=True)).fit()


# In[ ]:


print(thal_dmn_fit.summary())


# In[ ]:


print(thal_d_fit.summary())


# ## Hippocampal injections

# In[ ]:


hipp = structure_tree.get_structures_by_acronym(['HPF'])[0]
hipp_exps = mcc.get_experiments(cre=False, injection_structure_ids = [hipp['id']])
print(len(hipp_exps))


# In[ ]:


## Eliminate experiments with significant cortical leakage
leakage_experiments = []
for experiment in hipp_exps:
    ij_strs = [structure['id'] for structure in experiment['injection-structures']]
    ctx_strs_in_hipp_injection = [i for i in ij_strs if i in descendants]
    if len(ctx_strs_in_hipp_injection) > 0:
        leakage_experiments.append(experiment['id'])
print(len(leakage_experiments))


# In[ ]:


hipp_exps = [exp for exp in hipp_exps if exp['id'] not in leakage_experiments]
print(len(hipp_exps))


# In[ ]:


(hipp_centroids, hipp_injection_ratios_in, 
 hipp_projection_ratios_in, hipp_distances, 
 hipp_projections) = calculate_centroids_and_distances(hipp_exps, iso_mask)


# In[ ]:


hipp_d_coeff, hipp_dmn_coeff, hipp_tvals, hipp_pvals  = fit_glm(in_or_out_ctx, hipp_distances, hipp_projections)


# In[ ]:


hipp_exps = pd.DataFrame(hipp_exps)
hipp_glm_dat = pd.DataFrame({'id': hipp_exps['id'],
                           'injection structure': hipp_exps['structure-abbrev'],
                           'injection dmn fraction': hipp_injection_ratios_in,
                           'projection dmn fraction': hipp_projection_ratios_in,
                           'distance coefficient': hipp_d_coeff, 
                           'DMN coefficient': hipp_dmn_coeff,
                           't values': hipp_tvals,
                           'p values': hipp_pvals})


# In[ ]:


percent_dmn = []
dmn_coeff = []
for isid in hipp_exps['id']:
    percent_dmn.append(hipp_glm_dat[hipp_glm_dat['id'] == isid]['projection dmn fraction'].values[0] * 100)
    dmn_coeff.append(hipp_glm_dat[hipp_glm_dat['id'] == isid]['DMN coefficient'].values[0])
hipp_exps['percent_dmn'] = percent_dmn
hipp_exps['DMN_coefficient'] = dmn_coeff
hipp_exps = separate_coordinates(hipp_exps)
hipp_exps = get_injection_volume(hipp_exps)


# In[ ]:


hipp_glm_dat.to_csv(os.path.join(path, 'hippocampal injections distance dmn coefficients.csv'))


# In[ ]:


strs = mcc.rank_structures(hipp_glm_dat['id'], is_injection = True)
strids = [exp[0]['structure_id'] for exp in strs]
sorted_names = pd.DataFrame(
    hipp_glm_dat[['DMN coefficient', 
        'injection structure', 
        'injection dmn fraction', 
        'projection dmn fraction']].groupby(
        'injection structure').mean().sort_values(
        by = 'DMN coefficient', ascending = False))


# In[ ]:


fig, ax = plt.subplots(1, figsize = (4, 5))
sns.boxplot('injection structure', 'DMN coefficient', data = hipp_glm_dat, 
            order = sorted_names.index, ax = ax, palette = 'inferno_r')
ax.set_xticks(np.linspace(ax.get_xbound()[0]+0.5, ax.get_xbound()[1]-0.5, len(sorted_names.index)))
plt.xticks(rotation = -90)
plt.title('Hipppocampus Injections', y=1.03)
ax.axhline(y=0, xmin=0, xmax=100, color = 'grey', linestyle = 'dashed', zorder = -1)
ax.set_ylim(-0.5, 1)
plt.savefig(os.path.join(path, 'DMN_proj_by_source_wt_hippocampus_boxplot_horizontal.pdf'), 
            bbox_inches='tight', pad_inches=0.3, format='pdf', transparent=True, dpi=1000)


# ## Cortico-thalamic projections

# In[ ]:


ctx_exps = mcc.get_experiments(cre=False, injection_structure_ids = [iso['id']])
print(len(ctx_exps))


# In[ ]:


(ct_centroids, ct_injection_ratios_in, 
 ct_projection_ratios_in, ct_distances, 
 ct_projections) = calculate_centroids_and_distances(ctx_exps, th_mask)


# In[ ]:


in_or_out_th = shrink_to_mask(dmn_mask, th_mask)


# In[ ]:


ct_d_coeff, ct_dmn_coeff, ct_tvals, ct_pvals  = fit_glm(in_or_out_th, ct_distances, ct_projections)


# In[ ]:


ctx_exps = pd.DataFrame(ctx_exps)
ct_glm_dat = pd.DataFrame({'id': ctx_exps['id'],
                           'injection structure': ctx_exps['structure-abbrev'],
                           'injection dmn fraction': ct_injection_ratios_in,
                           'projection dmn fraction': ct_projection_ratios_in,
                           'distance coefficient': ct_d_coeff, 
                           'DMN coefficient': ct_dmn_coeff,
                           't values': ct_tvals,
                           'p values': ct_pvals})


# In[ ]:


percent_dmn = []
dmn_coeff = []
for isid in ctx_exps['id']:
    percent_dmn.append(ct_glm_dat[ct_glm_dat['id'] == isid]['projection dmn fraction'].values[0] * 100)
    dmn_coeff.append(ct_glm_dat[ct_glm_dat['id'] == isid]['DMN coefficient'].values[0])
ctx_exps['percent_dmn'] = percent_dmn
ctx_exps['DMN_coefficient'] = dmn_coeff
ctx_exps = separate_coordinates(ctx_exps)
ctx_exps = get_injection_volume(ctx_exps)


# In[ ]:


ct_glm_dat.to_csv(os.path.join(path, 'cortical injections distance dmn coefficients to thalamus.csv'))


# In[ ]:


strs = mcc.rank_structures(ct_glm_dat['id'], is_injection = True)
strids = [exp[0]['structure_id'] for exp in strs]
sorted_names = pd.DataFrame(
    ct_glm_dat[['DMN coefficient', 
        'injection structure', 
        'injection dmn fraction', 
        'projection dmn fraction']].groupby(
        'injection structure').mean().sort_values(
        by = 'DMN coefficient', ascending = False))


# In[ ]:


fig, ax = plt.subplots(1, figsize = (16, 6))
sns.boxplot('injection structure', 'DMN coefficient', data = ct_glm_dat, 
            order = sorted_names.index, ax = ax, palette = 'inferno_r')
ax.set_xticks(np.linspace(ax.get_xbound()[0]+0.5, ax.get_xbound()[1]-0.5, len(sorted_names.index)))
plt.xticks(rotation = -90)
plt.title('Cortical Injections', y=1.03)
ax.axhline(y=0, xmin=0, xmax=100, color = 'grey', linestyle = 'dashed', zorder = -1)
plt.savefig(os.path.join(path, 'DMN_proj_by_source_wt_cortical_injections_thalamus_projections_boxplot_horizontal.pdf'), 
            bbox_inches='tight', pad_inches=0.3, format='pdf', transparent=True, dpi=1000)


# ## plot injection locations for all wt injections: cortical, hippocampal, thalamic

# In[ ]:


mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                            resolution=25) # Switch to 25 um for images


# In[ ]:


all_exp_ids = list(ctx_exps['id'].values) + list(thal_exps['id'].values) + list(hipp_exps['id'].values)
len(all_exp_ids)


# In[ ]:


all_exps = pd.DataFrame(mcc.get_experiments(
    cre=False, 
    injection_structure_ids = [iso['id'], thal['id'], hipp['id']])) #25 um
all_exps = all_exps[all_exps['id'].isin(all_exp_ids)]
all_glm_dat = pd.concat([ctx_glm_dat, th_glm_dat, hipp_glm_dat])


# In[ ]:


percent_dmn = []
dmn_coeff = []
for isid in all_exps['id']:
    percent_dmn.append(all_glm_dat[all_glm_dat['id'] == isid]['projection dmn fraction'].values[0] * 100)
    dmn_coeff.append(all_glm_dat[all_glm_dat['id'] == isid]['DMN coefficient'].values[0])
all_exps['percent_dmn'] = percent_dmn
all_exps['DMN_coefficient'] = dmn_coeff


# In[ ]:


all_exps = separate_coordinates(all_exps)
all_exps = get_injection_volume(all_exps)


# In[ ]:


fig,ax = plt.subplots(1,1)
mask = get_mask('sagittal')
ax.imshow(template_sag, cmap = 'gray')
ax.imshow(mask, cmap = 'gist_gray', vmax = 3, alpha = 0.8)
ax.imshow(template_sag, cmap = 'gray', alpha = 0.6)
cax = ax.scatter(all_exps.zcoord/25,
           all_exps.ycoord/25,
           c = all_exps['DMN_coefficient'],
                cmap = 'inferno',
           edgecolor = 'none',
           s = all_exps.volume.values*500,
          alpha = 0.8)
cbar = fig.colorbar(cax, ax = ax, orientation='vertical', 
                    fraction=0.046, pad=0.04)
cbar_ax = fig.axes[-1]
cbar.solids.set_rasterized(True)
cbar.solids.set_edgecolor("face")
cbar.set_alpha(1)
cbar.draw_all()
ax.set_aspect('equal')
plt.axis('off')
fig.set_size_inches(10, 10)
plt.savefig(os.path.join(path, 'wt_all_inj_dmn_proj_injection_spatial_map_sag_DMN_coefficient_cmap_inferno.pdf'), 
           bbox_inches='tight', pad_inches=0.3, format='pdf', transparent = True, dpi=1000)


# ## Plot injection locations for cortical and hippocampal injections only

# In[ ]:


t_h_exps = pd.DataFrame(mcc.get_experiments(
    cre=False, 
    injection_structure_ids = [thal['id'], hipp['id']])) #25 um
good_exp_ids = list(thal_exps['id'].values) + list(hipp_exps['id'].values)
print(len(good_exp_ids))
t_h_exps = t_h_exps[t_h_exps['id'].isin(good_exp_ids)]
t_h_glm_dat = pd.concat([th_glm_dat, hipp_glm_dat])


# In[ ]:


percent_dmn = []
dmn_coeff = []
for isid in t_h_exps['id']:
    percent_dmn.append(t_h_glm_dat[t_h_glm_dat['id'] == isid]['projection dmn fraction'].values[0] * 100)
    dmn_coeff.append(t_h_glm_dat[t_h_glm_dat['id'] == isid]['DMN coefficient'].values[0])
t_h_exps['percent_dmn'] = percent_dmn
t_h_exps['DMN_coefficient'] = dmn_coeff


# In[ ]:


t_h_exps = separate_coordinates(t_h_exps)
t_h_exps = get_injection_volume(t_h_exps)


# In[ ]:


fig,ax = plt.subplots(1,1)
mask = get_mask('sagittal')
ax.imshow(template_sag, cmap = 'gray')
ax.imshow(mask, cmap = 'gist_gray', vmax = 3, alpha = 0.8)
ax.imshow(template_sag, cmap = 'gray', alpha = 0.6)
levels = np.arange(-0.75, 0.25, 1.25)
cax = ax.scatter(t_h_exps.zcoord/25,
           t_h_exps.ycoord/25,
           c = t_h_exps['DMN_coefficient'],
                cmap = 'inferno',
           edgecolor = 'none',
           s = t_h_exps.volume.values*500,
          alpha = 0.8,
                vmin = -0.75,
                vmax = 1.25)
cbar = fig.colorbar(cax, ax = ax, orientation='vertical', 
                    fraction=0.046, pad=0.04)
cbar_ax = fig.axes[-1]
cbar.solids.set_rasterized(True)
cbar.solids.set_edgecolor("face")
cbar.set_alpha(1)
cbar.draw_all()
ax.set_aspect('equal')
plt.axis('off')
fig.set_size_inches(10, 10)
plt.savefig(os.path.join(path, 'wt_hipp_thal_dmn_proj_injection_spatial_map_sag_DMN_coefficient_cmap_inferno.pdf'), 
           bbox_inches='tight', pad_inches=0.3, format='pdf', transparent = True, dpi=1000)


# ## Figure 3d

# In[117]:


print(ctx_d_fit.pvalues)
print(ctx_dmn_fit.pvalues)


# In[121]:


pvals = [ctx_dmn_fit.pvalues[1],
         ctx_d_fit.pvalues[1]]


# In[122]:


fdrcorr = sm.stats.fdrcorrection(pvals, alpha=0.05, method='indep')


# In[123]:


fdrcorr[1]


# In[124]:


ctx_glm_dat.keys()


# In[148]:


fig, ax = plt.subplots(figsize = (7, 7))
sns.regplot(x = 'injection dmn fraction',
            y = 'projection dmn fraction',
            data = ctx_glm_dat,
            ax = ax,
            color = 'purple',
            scatter_kws={'s':150})

ax.set_ylabel('Fraction of Cortical\nProjections in DMN')
ax.set_xlabel('Injection DMN Fraction')
ax.axhline(y = 0, xmin=0, xmax=100, color = 'k', linestyle = 'dashed')
plt.legend(loc = 2)
sns.despine()
ax.tick_params(top = False, right = False)

#everything on zorder -1 or lower will be rasterized
ax.set_rasterization_zorder(0)

plt.savefig(os.path.join(path, 'wt_injection_projection_DMN_fraction.pdf'), bbox_inches='tight', transparent=True, dpi=1000)


# In[150]:


fig, ax = plt.subplots(figsize = (7, 7))
ctx_glm_dat['scaled distance coefficient'] = ctx_glm_dat['distance coefficient'] * 10
ax1 = sns.regplot(x='injection dmn fraction', 
                  y='scaled distance coefficient', 
                  data = ctx_glm_dat,
                  ax = ax,
                  color = 'k',
                  label = 'Distance Coefficient, $R^2$={0} (p={1})'.format(round(ctx_d_fit.rsquared, 2),
                                                                           round(fdrcorr[1][1], 2)),
                  scatter_kws={'s':150})
sns.regplot(x = 'injection dmn fraction',
            y = 'DMN coefficient',
            data = ctx_glm_dat,
            ax = ax,
            color = 'orange',
            label = 'DMN coefficient, $R^2$={0} (p={1:0.1e})'.format(round(ctx_dmn_fit.rsquared, 2),
                                                             fdrcorr[1][0]),
            scatter_kws={'s':150})

ax.set_ylabel('DMN Coefficient,\nDistance Coefficient x 10')
ax.set_xlabel('Injection DMN Fraction')
ax.axhline(y = 0, xmin=0, xmax=100, color = 'k', linestyle = 'dashed')
ax.set_xlim([-.2, 1.2])
plt.legend(loc = 2, fontsize = 15)
sns.despine()
ax.tick_params(top = False, right = False)

#everything on zorder -1 or lower will be rasterized
ax.set_rasterization_zorder(0)

plt.savefig(os.path.join(path, 'wt_distance_DMN_coeff.pdf'), bbox_inches='tight', transparent=True, dpi=1000)


# In[ ]:


fig, ax = plt.subplots(figsize = (7, 7))
th_glm_dat['scaled distance coefficient'] = th_glm_dat['distance coefficient'] * 10
ax1 = sns.regplot(x='injection dmn fraction', 
                  y='scaled distance coefficient', 
                  data = th_glm_dat,
                  ax = ax,
                  color = 'k',
                  label = 'Distance Coefficient, $R^2$={0} (p={1})'.format(round(thal_d_fit.rsquared, 2),
                                                                           round(fdrcorr[1][3], 1)),
                  scatter_kws={'s':150})
sns.regplot(x = 'injection dmn fraction',
            y = 'DMN coefficient',
            data = th_glm_dat,
            ax = ax,
            color = 'orangered',
            label = 'DMN coefficient, $R^2$={0} (p={1})'.format(round(thal_dmn_fit.rsquared, 3),
                                                             round(fdrcorr[1][2], 1)),
            scatter_kws={'s':150})

ax.set_ylabel('DMN Coefficient, Distance Coefficient x 10')
ax.set_xlabel('Injection DMN Fraction')
ax.axhline(y = 0, xmin=0, xmax=100, color = 'k', linestyle = 'dashed')
ax.set_ylim([-2, 1.5])
plt.legend(loc = 2)

#everything on zorder -1 or lower will be rasterized
ax.set_rasterization_zorder(0)

plt.savefig(os.path.join(path, 'wt_distance_DMN_coeff_thalamus_injections.pdf'), bbox_inches='tight', transparent=True, dpi=1000)


# In[ ]:


ctx_glm_dat.head()


# In[ ]:


print(len(ctx_glm_dat[ctx_glm_dat['injection dmn fraction'] == 1]))
print(1-ctx_glm_dat[ctx_glm_dat['injection dmn fraction'] == 1]['projection dmn fraction'].min())
print(1-ctx_glm_dat[ctx_glm_dat['injection dmn fraction'] == 1]['projection dmn fraction'].max())
print(len(ctx_glm_dat[ctx_glm_dat['injection dmn fraction'] > 0.95]))
print(1-ctx_glm_dat[ctx_glm_dat['injection dmn fraction'] > 0.95]['projection dmn fraction'].min())
print(1-ctx_glm_dat[ctx_glm_dat['injection dmn fraction'] > 0.95]['projection dmn fraction'].max())


# In[ ]:




