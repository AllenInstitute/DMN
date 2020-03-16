# -*- coding: utf-8 -*-
"""
Spyder Editor

Enter the image series id of interest (target-defined) to plot the layer
distribution for that experiment and all matched controls
"""
import os
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_context('paper')
sns.set_style('white')

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import platform
if platform.system() == 'Darwin':
    path = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/_new_figures/Figure_5'
    dat = pd.read_csv(r'/Users/jenniferwh/Dropbox/DMN data/correlations/_final/good_td_wt_correlations_with_inj_corr.csv')
    metadat = pd.read_csv(r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/target_defined_dataset.csv',
                          engine = 'python')
elif platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\_new_figures\Figure_5' 
    dat = pd.read_csv(r'C:\Users\jenniferwh\Dropbox (Personal)\DMN data\correlations\_final\good_td_wt_correlations_with_inj_corr.csv')
    metadat = pd.read_csv(r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\target_defined_dataset.csv',
                          engine = 'python')
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file = 'connectivity/mouse_connectivity_manifest.json')

st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
ai_map = {value:key for key,value in ia_map.items()}
st_dict = dict(zip(metadat['image_series_id'], + metadat['source']+ '$_{' + 
                   metadat['target_by_projection'] + '}$'))
structs = st.get_structures_by_set_id(
        [667481440, 667481441, 667481445, 667481446, 667481449]
        )
acronyms = [structure['acronym'] for structure in structs]
ids = [structure['id'] for structure in structs]
ids += [1121, 526, 20, 52, 543, 664, 92, 712, 139, 727, 28, 60, 743] #ENTl, ENTm
L1strs = [ia_map[structure] for structure in acronyms if '1' in structure]
L1strs += [526, 1121]
L2strs = [ia_map[structure] for structure in acronyms if '2/3' in structure]
L2strs += [20, 543, 52, 664]
L4strs = [ia_map[structure] for structure in acronyms if '4' in structure]
L4strs += [92, 712]
L5strs = [ia_map[structure] for structure in acronyms if '5' in structure]
L5strs += [139, 727]
L6strs = [ia_map[structure] for structure in acronyms if '6' in structure]
L6strs = [structure for structure in L6strs if '6b' not in ai_map[structure]] #6b is unreliable
L6strs += [28, 743]

# remove L6b structures from quantified injection area
L6bstrs = [structure for structure in L6strs if '6b' in ai_map[structure]]
structs = [structure for structure in structs if structure not in L6bstrs]

def plot_clustered_stacked(dfall, labels=None, title = None, H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe
https://stackoverflow.com/questions/22787209/how-to-have-clusters-of-stacked-bars-with-python-pandas"""
    fig = plt.figure()
    pad = 16
    n_df = len(dfall)
    n_col = len(dfall[0].columns)
    n_ind = len(dfall[0].index)
    axe = plt.subplot(111)
    layer_colors = ['k', 'c', 'm', 'gray', 'gold']
    
    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() / float(n_df + pad) * i )
                rect.set_width(1 / float(n_df + pad))
                rect.set_color(layer_colors[j])
    minval = 0.02
    maxval = 0.22
    step = (maxval-minval)/n_df
    axe.set_xticks(np.arange(minval, maxval, step))
    axe.tick_params(axis='both', which='major', pad=-3)
    axe.set_ylim([0.01, 1]) #this gets rid of horizontal line on top of plot
    axe.set_yticklabels([None])
    axe.set_xlabel("")
    axe.set_xticklabels(labels, rotation = -30, fontsize = 8)
    axe.set_title(title, fontsize = 16,
                  fontweight='bold')
    axe.invert_yaxis()
    
    axe.spines['top'].set_visible(False)
    axe.spines['right'].set_visible(False)
    axe.spines['bottom'].set_visible(False)
    axe.spines['left'].set_visible(False)

    # Add invisible data to add another legend
    n=[]   
     
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=1, fontsize = 6,
                    bbox_to_anchor = [0.8,1],
                    labelspacing = 0.2)
    axe.add_artist(l1)
    fig.set_size_inches(2.5, 1.7)
    plt.tight_layout()
    return axe

isids = [609475867, 475829896, 649362978, 637855050, 477435412,
         571401645, 528741104, 607321130, 607059419] #right to left in figure
df = pd.DataFrame({'experiment': isids})
injs = mcc.get_structure_unionizes(experiment_ids = df['experiment'].values, 
                                       is_injection = True, hemisphere_ids = [3]) 
for ix, row in df.iterrows():
    print(row['experiment'])
    inj_total = injs[(injs['experiment_id'] == row['experiment']) & 
                 (injs['structure_id'].isin(ids))
                 ]['projection_volume'].sum()
    df.loc[df['experiment'] == row['experiment'], 
            'L1'] = injs[(injs['experiment_id'] == row['experiment']) & 
                         (injs['structure_id'].isin(L1strs))
                         ]['projection_volume'].sum()/inj_total
    df.loc[df['experiment'] == row['experiment'], 
            'L2/3'] = injs[(injs['experiment_id'] == row['experiment']) & 
                         (injs['structure_id'].isin(L2strs))
                         ]['projection_volume'].sum()/inj_total
    df.loc[df['experiment'] == row['experiment'], 
            'L4'] = injs[(injs['experiment_id'] == row['experiment']) & 
                         (injs['structure_id'].isin(L4strs))
                         ]['projection_volume'].sum()/inj_total
    df.loc[df['experiment'] == row['experiment'], 
            'L5'] = injs[(injs['experiment_id'] == row['experiment']) & 
                           (injs['structure_id'].isin(L5strs))
                           ]['projection_volume'].sum()/inj_total
    df.loc[df['experiment'] == row['experiment'], 
            'L6'] = injs[(injs['experiment_id'] == row['experiment']) & 
                       (injs['structure_id'].isin(L6strs))
                       ]['projection_volume'].sum()/inj_total
dfs = []
for ix, row in df.iterrows():
    dfs.append(pd.DataFrame(row).transpose().set_index('experiment'))
#dfs= sorted(dfs,key=lambda x:x["L6"].max(axis=0))
#isids = [int(df.index.values[0]) for df in dfs]
try:
    labels = [st_dict[isid] for isid in isids]
except:
    labels = isids
fname = 'ACAd_all_TD'
labels = ['ACAd$_{RSPagl}$', 'ACAd$_{RSPv}$', 'ACAd$_{RSPv}$', 'ACAd$_{RSPv}$',
          'ACAd$_{VISam}$', 'ACAd$_{VISam}$', 'ACAd$_{ORBvl}$', 'ACAd$_{MOp}$',
          'ACAd$_{VISp}$'] #left to right in figure
g = plot_clustered_stacked(dfs, labels)

plt.savefig(os.path.join(path, 
                         'layer_quantification_{0}.pdf'.format(fname)),
            bbox_inches='tight', 
            pad_inches=0.3, format='pdf', transparent = True, dpi=300)

    