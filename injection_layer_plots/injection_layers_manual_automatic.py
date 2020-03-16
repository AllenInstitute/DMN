# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_context('paper')
sns.set_style('white')

import matplotlib.pyplot as plt

import platform
if platform.system() == 'Darwin':
    path = r'/Users/jenniferwh/Dropbox/DMN data/layers'
elif platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox (Personal)\DMN data\layers' 
    
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file = 'connectivity/mouse_connectivity_manifest.json')

dat = pd.read_csv(os.path.join(path, 'manual_layer_counts_sig_experiments.csv'))
dat.loc[dat['L6'].isnull(), 'L6'] = 0
st_dict = dict(zip(dat['experiment'], dat['source']+ '/' + dat['target']))
injs = mcc.get_structure_unionizes(experiment_ids = dat['experiment'].values, is_injection = True)

st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
ai_map = {value:key for key,value in ia_map.items()}
structs = st.get_structures_by_set_id(
        [667481440, 667481441, 667481445, 667481446, 667481449, 667481450])
acronyms = [structure['acronym'] for structure in structs]
ids = [structure['id'] for structure in structs]
ids += [1121, 526, 20, 52, 543, 664, 92, 712, 139, 727, 28, 60, 743] #ENTl, ENTm
L1strs = [ia_map[structure] for structure in acronyms if '1' in structure]
L1strs += [526, 1121]
L2strs = [ia_map[structure] for structure in acronyms if '2/3' in structure]
L2strs += [20, 543]
L3strs = [52, 664]
L5strs = [ia_map[structure] for structure in acronyms if '5' in structure]
L5strs += [139, 727]
L6strs = [ia_map[structure] for structure in acronyms if '6' in structure]
L6strs = [structure for structure in L6strs if '6b' not in ai_map[structure]] #6b is unreliable
L6strs += [28, 743]

# remove L6b structures from quantified injection area
L6bstrs = [structure for structure in L6strs if '6b' in ai_map[structure]]
structs = [structure for structure in structs if structure not in L6bstrs]

for row in dat.iterrows():
    total_cells = np.sum(row[1][['L2/3', 'L5', 'L6']])
    dat.loc[dat['experiment'] == row[1]['experiment'], 'L2/3 ratio'] = row[1]['L2/3']/total_cells
    dat.loc[dat['experiment'] == row[1]['experiment'], 'L5 ratio'] = row[1]['L5']/total_cells
    dat.loc[dat['experiment'] == row[1]['experiment'], 'L6 ratio'] = row[1]['L6']/total_cells
    
    relevant_structs = L2strs + L5strs + L6strs
    inj_total = injs[(injs['experiment_id'] == row[1]['experiment']) &
                     (injs['structure_id'].isin(relevant_structs))]['projection_volume'].sum()
    dat.loc[dat['experiment'] == row[1]['experiment'], 
            'L2/3 actual'] = injs[(injs['experiment_id'] == row[1]['experiment']) &
                         (injs['structure_id'].isin(L2strs))]['projection_volume'].sum()/inj_total
    dat.loc[dat['experiment'] == row[1]['experiment'], 
            'L5 actual'] = injs[(injs['experiment_id'] == row[1]['experiment']) &
                         (injs['structure_id'].isin(L5strs))]['projection_volume'].sum()/inj_total
    dat.loc[dat['experiment'] == row[1]['experiment'], 
            'L6 actual'] = injs[(injs['experiment_id'] == row[1]['experiment']) &
                         (injs['structure_id'].isin(L6strs))]['projection_volume'].sum()/inj_total

def plot_clustered_stacked(dfall, labels=None, title="multiple stacked bar plot",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe
https://stackoverflow.com/questions/22787209/how-to-have-clusters-of-stacked-bars-with-python-pandas"""
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colormap = [colors[4], colors[3], colors[1]]
    
    n_df = len(dfall)
    n_col = len(dfall[0].columns) -1
    n_ind = len(dfall[0].index)
    axe = plt.subplot(111)
    for ix, df in enumerate(dfall) : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      color = colormap,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.label, rotation = 90) #df.index for isid, df.label for text
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]        
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    axe.add_artist(l1)
    return axe


meandata = dat[['experiment', 'L6 ratio', 'L5 ratio', 'L2/3 ratio']].melt(id_vars = 'experiment', 
              var_name = 'layer')
meandata['method'] = 'manual'
meandatb = dat[['experiment', 'L6 actual',  'L5 actual', 'L2/3 actual']].melt(id_vars = 'experiment',
              var_name = 'layer')
meandatb['method'] = 'automatic'
meandat = pd.concat([meandata, meandatb])
for layer in [ 'L6', 'L5', 'L2/3']:
    meandat.loc[meandat['layer'].str.contains(layer), 'layer'] = layer
pivot_dfm = meandat[meandat['method'] == 'manual'].pivot(
        index = 'experiment', 
        columns = 'layer', 
        values = 'value')
pivot_dfm['label'] = [st_dict[experiment] for experiment in pivot_dfm.index]
pivot_dfm.sort_values(by='label', inplace = True)
pivot_dfm = pivot_dfm[['L6', 'L5', 'L2/3', 'label']]
pivot_dfa = meandat[meandat['method'] == 'automatic'].pivot(
        index = 'experiment', 
        columns = 'layer', 
        values = 'value')
pivot_dfa['label'] = [st_dict[experiment] for experiment in pivot_dfa.index]
pivot_dfa.sort_values(by='label', inplace = True)
pivot_dfa = pivot_dfa[['L6', 'L5', 'L2/3', 'label']]
g = plot_clustered_stacked([pivot_dfm, pivot_dfa],["manual", "automatic"])

g.axes.yaxis.set_ticklabels([])
plt.yticks([])   
plt.savefig(os.path.join(path, 'layer_quantification_manual_automated.png'),
            bbox_inches='tight', 
            pad_inches=0.3, format='png', transparent = True, dpi=300)
    