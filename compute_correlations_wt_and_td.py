# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
import os
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi
import scipy.stats as stats
from scipy.optimize import curve_fit
import platform
if platform.system() == 'Darwin':
    path = '/Users/jenniferwh/Dropbox/DMN data/correlations'
    outpath = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN'
elif platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox\DMN data\correlations'
    outpath = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN'
savepath = os.path.join(outpath, 'Figure 7 td-wt correlations by source')

#%% Data curation
# Start wtih all correlations: output of cluster code, but manually changed 
# "Fraction of match for <60% primary" in Excel

aapi = AnatomyApi()
ss = aapi.get_summary_structure_data('id')
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
structure_tree = mcc.get_structure_tree()
isocortex = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
cla = structure_tree.get_structures_by_acronym(['CLA'])[0]['id']
HPF = structure_tree.get_structures_by_acronym(['HPF'])[0]
iso = structure_tree.descendant_ids([isocortex['id']])[0]
iso = [structure for structure in iso if structure in ss]
hipp = structure_tree.descendant_ids([HPF['id']])[0]
hipp = [structure for structure in hipp if structure in ss]
ia_map = structure_tree.get_id_acronym_map()
ai_map = {value:key for key, value in ia_map.items()}
ctx_strs = [ai_map[structure] for structure in iso]
hipp_strs = [ai_map[structure] for structure in hipp]
valid_strs = ctx_strs+hipp_strs+[ai_map[cla]]

td_dataset = pd.read_csv(os.path.join(outpath, 'target_defined_dataset.csv'))
td_dataset = td_dataset[td_dataset['include'] == 'yes']
td_dataset = td_dataset[td_dataset['source'].isin(valid_strs)]

distance_threshold = 800
overlap_threshold = 0.01
alldat = pd.read_csv(os.path.join(path, 'td_wt_cre_matched_correlations_NPV_ipsi_with_inj_corr.csv'))
alldat = alldat[alldat['image_series_id'].isin(td_dataset['image_series_id'])]
c_by_source = pd.read_csv(os.path.join(path, 'match_correlations_by_source_NPV_ipsi_with_inj_corr.csv'))
c_by_source = c_by_source[c_by_source['same_primary'] == True]
c_by_source = c_by_source[c_by_source['same secondary for <60% primary'] != False]
alldat = alldat[alldat['same_primary'] == True]
alldat = alldat[alldat['same secondary for <60% primary'] != False]
# remove duplicates
c_by_source['index_original'] = c_by_source.groupby(['match_A', 'match_B']).match_A.transform('idxmin')    
c_by_source = c_by_source[~c_by_source.duplicated(subset=['match_A', 'match_B'], keep='first')]
for isid in c_by_source['match_A'].unique():
    Bmatches = c_by_source[c_by_source['match_A'] == isid]['match_B'].values
    Amatches = c_by_source[c_by_source['match_B'] == isid]['match_A'].values
    duplicates = [match for match in Amatches if match in Bmatches]
    if len(duplicates) > 0:
        print(Bmatches)
        print(Amatches)

alldat = alldat.merge(td_dataset[['image_series_id', 'include']], 
                      on='image_series_id', 
                      how = 'left')
alldat = alldat[alldat['include'] == 'yes']

# wild type experiments excluded
# experiments with leakage, tile edges, low GFP signal, surface artifacts, 
# wrong injection site assigned, not all layers labeled, 
# or looks like PT cells only
fail_expts = [114008926, 120280939, 180073473, 180403712, 180601025, 183174303, 183329222,
              249396394, 296047806, 299446445, 301060890, 303784745, 480069939, 482578964, 
              506947040, 514333422, 525796603, 545428296, 559878074, 638314843, 182888003,
             304585910, 183171679, 272930013, 523718075, 517072832, 148964212, 304762965,
             566992832, 272930013, 304762965, 266250904, 114399224, 286483411, 286417464,
             593277684, 546103149, 642809043, 286483411, 304564721] #VISp outlier excluded

alldat = alldat[~alldat['match_id'].isin(fail_expts)]
c_by_source = c_by_source[~c_by_source['match_A'].isin(fail_expts)]
c_by_source = c_by_source[~c_by_source['match_B'].isin(fail_expts)]
alldat = alldat[alldat['distance'] < distance_threshold]
c_by_source = c_by_source[c_by_source['distance'] < distance_threshold]

c_by_source = c_by_source[c_by_source['fraction of match covered by td injection'] > overlap_threshold]
alldat = alldat[(alldat['fraction of match covered by td injection'] > overlap_threshold) |
               (alldat['fraction of match covered by td exclusion zone'] > overlap_threshold)]
alldat_extra = pd.read_csv(os.path.join(path, 'td_wt_cre_matched_correlations_NPV_ipsi_with_inj_corr.csv'))
alldat_extra = alldat_extra[alldat_extra['JW_pass-fail'] == 'Include as exception']
alldat = pd.concat([alldat, alldat_extra], sort=True)

# Remove thalamus and hippocampus proper sources for this analysis
c_by_source = c_by_source[(c_by_source['match_A_primary_source'].isin(valid_strs)) |
                    (c_by_source['match_B_primary_source'].isin(valid_strs))]
alldat = alldat[alldat['source'].isin(valid_strs)]
c_by_source['source'] = c_by_source['match_A_primary_source']
c_by_source = c_by_source[~c_by_source['source'].isin(['CA1', 'CA3', 'DG', 'SUB'])]
alldat = alldat[~alldat['source'].isin(['CA1', 'CA3', 'DG', 'SUB'])]

print(len(td_dataset))
print(len(alldat))

#%% model fit (starting at best model)
def exp(X, a, b, c, d, e):
    x,y,z = X
    return a * 10**(-x / b) + (c * y) + (d * z) + e
c_by_source.sort_values(by='match_A_injection_size', inplace = True)
x = c_by_source['match_A_injection_size'].values
y = c_by_source['distance'].values
y2 = c_by_source['fraction of match covered by td injection'].values
Y = c_by_source['spearman_correlation'].values

popt, pcov = curve_fit(exp, (x,y,y2), Y)

yfit=exp((x,y,y2), *popt)

resids = yfit-Y

ssq1=((yfit-Y)**2).sum()

n = len(x)    # number of data points
p = len(popt) # number of parameters

df = max(0, n-p) # number of degrees of freedom

ss_res = np.sum((Y - yfit) ** 2)

# total sum of squares
ss_tot = np.sum((Y - np.mean(Y)) ** 2)

# r-squared
r2 = 1 - (ss_res / ss_tot)
print(r2)
print(ss_res)
print(popt)

c_by_source['exp_predicted'] = exp((c_by_source['match_A_injection_size'].values,
                               c_by_source['distance'].values,
                                   c_by_source['fraction of match covered by td injection'].values), *popt)
alldat['exp_predicted'] = exp((alldat['td_injection_size'],
                               alldat['distance'].values,
                              alldat['fraction of match covered by td injection'].values), *popt)

#%% 95% prediction band
def predband(x, xd, yd, f_vars, conf=0.95):
    """
    Code adapted from Rodrigo Nemmen's post:
    http://astropython.blogspot.com.ar/2011/12/calculating-prediction-band-
    of-linear.html

    Calculates the prediction band of the regression model at the
    desired confidence level.

    Clarification of the difference between confidence and prediction bands:

    "The prediction bands are further from the best-fit line than the
    confidence bands, a lot further if you have many data points. The 95%
    prediction band is the area in which you expect 95% of all data points
    to fall. In contrast, the 95% confidence band is the area that has a
    95% chance of containing the true regression line."
    (from http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?
    reg_graphing_tips_linear_regressio.htm)

    Arguments:
    - x: array with x values to calculate the confidence band.
    - xd, yd: data arrays.
    - a, b, c: linear fit parameters.
    - conf: desired confidence level, by default 0.95 (2 sigma)

    References:
    1. http://www.JerryDallal.com/LHSP/slr.htm, Introduction to Simple Linear
    Regression, Gerard E. Dallal, Ph.D.
    """
    for ix, val in enumerate(xd):
        if all(x == val):
            index = ix
    alpha = 1. - conf    # Significance
    N = len(xd) * len(xd[index])          # data sample size
    var_n = len(f_vars)  # Number of variables used by the fitted function.

    # Quantile of Student's t distribution for p=(1 - alpha/2)
    q = stats.t.ppf(1. - alpha / 2., N - var_n)

    # Std. deviation of an individual measurement (Bevington, eq. 6.15)
    se = np.sqrt(1. / (N - var_n) * np.sum((yd - exp(xd, *f_vars)) ** 2))

    # Auxiliary definitions
    sx = (x - xd[0].mean()) ** 2
    sxd = np.sum((xd[0] - xd[index].mean()) ** 2)

    # Predicted values (best-fit model)
    yp = exp(xd, *f_vars)
    # Prediction band
    dy = q * se * np.sqrt(1. + (1. / N) + (sx / sxd))

    # Upper & lower prediction bands.
    lpb, upb = yp - dy, yp + dy

    return lpb, upb

low_pred = pd.DataFrame()
for key in ['td_injection_size', 'distance', 'fraction of match covered by td injection']:
    low_pred[key] = predband(alldat[key].values, 
                                       (alldat['td_injection_size'].values, 
                                        alldat['distance'].values,
                                        alldat['fraction of match covered by td injection'].values), 
                                       alldat['spearman_correlation'].values, popt)[0]
                                       
alldat['low_pred_band'] = low_pred.min(axis = 1).values
alldat['sig'] = alldat['spearman_correlation'] < alldat['low_pred_band']

low_pred = pd.DataFrame()
for key in ['match_A_injection_size', 'distance', 'fraction of match covered by td injection']:
    low_pred[key] = predband(c_by_source[key].values, 
                                       (c_by_source['match_A_injection_size'].values, 
                                        c_by_source['distance'].values,
                                        c_by_source['fraction of match covered by td injection'].values), 
                                       c_by_source['spearman_correlation'].values, popt)[0]
c_by_source['low_pred_band'] = low_pred.min(axis = 1).values
c_by_source['sig'] = c_by_source['spearman_correlation'] < c_by_source['low_pred_band']
c_by_source['experiment'] = 'Control'
alldat['experiment'] = 'Target-Defined'

#%% add DMN and contralateral correlation data
core_dmn = ['ACAv', 'ACAd', 'SSp-tr', 'SSp-ll', 'ILA', 'PL', 'ORBvl', 'VISa', 'ORBm', 'ORBl', 'VISam', 'RSPv',
            'RSPd', 'MOs', 'RSPagl']
dmn = core_dmn + ['VISrl', 'VISpm', 'AId', 'AIv']
out_dmn = ['VISp', 'VISpor', 'VISal', 'SSp-n', 'AUDpo', 'SSp-bfd', 'VISl', 'AUDp', 'AUDd', 'SSp-m', 'MOp', 'SSs',
           'SSp-ul', 'VISC', 'ECT', 'VISli', 'TEa']
c_by_source.loc[c_by_source['source'].isin(['RSPagl', 'VISrl', 'VISpm', 'AIv']), 'group'] = 'DMN'
c_by_source.loc[c_by_source['source'].isin(out_dmn), 'group'] = 'out-DMN'
c_by_source.loc[c_by_source['source'].isin(core_dmn), 'group'] = 'DMN Core'
alldat.loc[alldat['source'].isin(['RSPagl', 'VISrl', 'VISpm', 'AIv']), 'group'] = 'DMN'
alldat.loc[alldat['source'].isin(out_dmn), 'group'] = 'out-DMN'
alldat.loc[alldat['source'].isin(core_dmn), 'group'] = 'DMN Core'

# ic = ipsi contra
alldat_ic = pd.read_csv(os.path.join(path, r'matched_td_data_with_sig.csv'))

alldat_ic.rename(columns = {'spearman_correlation': 'spearman_correlation_with_contra', 
                           'sig': 'sig_with_contra'}, inplace = True)
alldat = alldat.merge(alldat_ic[['image_series_id', 'match_id', 
                                 'spearman_correlation_with_contra',
                                'sig_with_contra']], on = ['image_series_id', 'match_id'], how = 'left')

alldat.to_csv(os.path.join(path, 'good_td_wt_correlations_with_inj_corr.csv'), index = False)
c_by_source.to_csv(os.path.join(path, 'good_wt_correlations_with_inj_corr.csv'), index = False)
#%% Means
# Wild type matches
all_isids = np.unique(np.concatenate((c_by_source['match_A'].unique(), c_by_source['match_B'].unique())))
num_comparisons = []
corr = []
predicted_corr = []
sources = []
inj_sizes = []
mean_ol = []
mean_dist = []
mean_low_pred = []
for isid in all_isids:
    if len(c_by_source[c_by_source['match_A'] == isid]['match_A_primary_source'].values) > 0:
        sources.append(c_by_source[c_by_source['match_A'] == isid]['match_A_primary_source'].values[0])
        inj_sizes.append(c_by_source[c_by_source['match_A'] == isid]['match_A_injection_size'].values[0])
    else:
        sources.append(c_by_source[c_by_source['match_B'] == isid]['match_B_primary_source'].values[0])
        inj_sizes.append(c_by_source[c_by_source['match_B'] == isid]['match_A_injection_size'].values[0])
    dataset = c_by_source[(c_by_source['match_A'] == isid) | (c_by_source['match_B'] == isid)]
    num_comparisons.append(len(dataset))
    corr.append(dataset['spearman_correlation'].mean())
    predicted_corr.append(dataset['exp_predicted'].mean())
    mean_ol.append(dataset['injection_overlap'].mean())
    mean_dist.append(dataset['distance'].mean())
    mean_low_pred.append(dataset['low_pred_band'].min())
cmeans = pd.DataFrame({'image_series_id': all_isids, 'source': sources, 
                       'injection_size': inj_sizes, 'number_comparisons': num_comparisons, 
                       'spearman_correlation': corr, 'predicted_spearman': predicted_corr,
                      'mean_overlap': mean_ol, 'mean_distance': mean_dist, 'low_pred': mean_low_pred})
print(len(cmeans))
cmeans.to_csv(os.path.join(savepath, 'control_mean_corr_dat_ipsi.csv'),
              index = False)

#target-defined - wt matches
meandat = alldat.groupby('image_series_id')[['spearman_correlation', 
                                             'exp_predicted', 'low_pred_band', 'inj_corr']].mean().reset_index()
meandat['number of comparisons'] = alldat.groupby('image_series_id')['match_id'].count().values
for isid in meandat['image_series_id']:
    dataset = alldat[alldat['image_series_id'] == isid]
    meandat.loc[meandat['image_series_id'] == isid, 'source'] = alldat[
        alldat['image_series_id'] == isid]['source'].unique()[0]
    meandat.loc[meandat['image_series_id'] == isid, 'injection_size'] = alldat[
        alldat['image_series_id'] == isid]['td_injection_size'].unique()[0]
    meandat.loc[meandat['image_series_id'] == isid, 'mean_overlap'] = alldat[
        alldat['image_series_id'] == isid]['injection_overlap'].mean()
    meandat.loc[meandat['image_series_id'] == isid, 'low_pred_band'] = alldat[
        alldat['image_series_id'] == isid]['low_pred_band'].min()
    
meandat = meandat.merge(td_dataset[['image_series_id', 'target_by_projection']], 
                        on = 'image_series_id',
                        how = 'left')

meandat_ic = alldat_ic.groupby('image_series_id')[['spearman_correlation_with_contra', 
                                             'exp_predicted', 'low_pred_band']].mean().reset_index()
meandat_ic['number of comparisons'] = alldat_ic.groupby('image_series_id')['match_id'].count().values
for isid in meandat_ic['image_series_id']:
    dataset = alldat_ic[alldat_ic['image_series_id'] == isid]
    meandat_ic.loc[meandat_ic['image_series_id'] == isid, 'source'] = alldat_ic[
        alldat_ic['image_series_id'] == isid]['source'].unique()[0]
    meandat_ic.loc[meandat_ic['image_series_id'] == isid, 'injection_size'] = alldat_ic[
        alldat_ic['image_series_id'] == isid]['td_injection_size'].unique()[0]
    meandat_ic.loc[meandat_ic['image_series_id'] == isid, 'mean_overlap'] = alldat_ic[
        alldat_ic['image_series_id'] == isid]['injection_overlap'].mean()
meandat = meandat.merge(meandat_ic[['image_series_id', 'spearman_correlation_with_contra']],
                       on = 'image_series_id', how = 'left')

meandat.to_csv(os.path.join(savepath, 'td_mean_corr_dat_ipsi.csv'), 
               index = False)





