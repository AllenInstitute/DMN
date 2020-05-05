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
import json
sns.set_context('paper')
sns.set_style('white')

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import platform
if platform.system() == 'Darwin':
    path = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/_new_figures/Figure_5'
    dat = pd.read_csv(r'/Users/jenniferwh/Dropbox/DMN data/correlations/_final/good_td_wt_correlations_with_inj_corr.csv')
    metadat = pd.read_csv(r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/target_defined_dataset.csv')
elif platform.system() == 'Windows':
    path = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\_new_figures\Figure_5' 
    dat = pd.read_csv(r'C:\Users\jenniferwh\Dropbox (Personal)\DMN data\correlations\_final\good_td_wt_correlations_with_inj_corr.csv')
    metadat = pd.read_csv(r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN\target_defined_dataset.csv')
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache(manifest_file = 'connectivity/mouse_connectivity_manifest.json')

st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
ai_map = {value:key for key,value in ia_map.items()}
st_dict = dict(zip(metadat['image_series_id'], + metadat['source']+ '$_{' + 
                   metadat['target_by_polygon'] + '}$'))
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
    layer_colors = ['k', '#00ffff', '#ff00ff', '#008000', '#ffff00']
    
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
    minval = -0.38
    maxval = 0.08
    step = (maxval-minval)/n_df
    print(step)
    axe.set_xticks(np.arange(minval, maxval, step))
    #axe.set_xticks([-0.05, 0.02, 0.12])
    axe.tick_params(axis='both', which='major', pad=-3)
    #This may hide a little bit of the top of the bars - use with caution.
    #axe.set_ylim([0.01, 1]) #this gets rid of horizontal line on top of plot
    axe.set_yticklabels([None])
    axe.set_xlabel("")
    axe.set_xticklabels(labels, rotation = 90, fontsize = 6)
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

    l1 = axe.legend(h[:n_col], l[:n_col], loc=1, fontsize = 5,
                    bbox_to_anchor = [0.8,1],
                    labelspacing = 0.2)
    axe.add_artist(l1)
    fig.set_size_inches(1.75, 1.7)
    plt.tight_layout()
    return axe

all_isids = [112595376, 666090944, 569904687, 868641659, 592724077, 561511939, 521255975, #868641659
         623838656, 592522663] #right to left in figure
isids = [112595376, 666090944, 569904687, 592724077, 561511939, 521255975, #868641659
         623838656, 592522663] #available through SDK
not_online = [868641659]
df = pd.DataFrame({'experiment': all_isids})

inj_unionizes = []
if len(not_online)>0: # new data not released online
    for expt in not_online:
        if platform.system() == 'Windows':
            unionize_path = os.path.join(r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\DMN_paper\alternative_unionizes',
                                         'experiment_{0}'.format(str(expt)), 
                                         'output.json') #new data not online yet
        else:
            unionize_path = os.path.join(r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/data_files/alternative_unionizes',
                                         'experiment_{0}'.format(str(expt)), 
                                         'output.json') #new data not online yet
        with open(unionize_path, 'r') as jsonfile:
            unionize_dat = json.load(jsonfile)
        unionize_dat = pd.DataFrame(unionize_dat)
        unionize_dat = unionize_dat[unionize_dat['threshold'] == 0]
        unionize_dat.drop_duplicates(inplace = True)
        unionize_dat.rename(columns = {'projection_volume_above_threshold': 'projection_volume',
                                   'normalized_projection_volume_above_threshold': 
                                       'normalized_projection_volume',
                                       'image_series_id': 'experiment_id'}, 
            inplace = True)
        inj_unionize = unionize_dat[(unionize_dat['is_injection'] == True) &
                              (unionize_dat['hemisphere_id'] == 3)]
        inj_unionizes.append(inj_unionize)

injs = mcc.get_structure_unionizes(experiment_ids = isids, 
                                   is_injection = True, hemisphere_ids = [3]) 
injs = pd.concat([injs, inj_unionize], sort = True)
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
fname = 'RSPv'
labels = ['PL', 'PL', 'ACAd', 'ACAd', 'VISpl', 'VISp', 'VISl', 'VISpl', 'WT'] #left to right in figure
g = plot_clustered_stacked(dfs, labels)

savepath = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN/_new_figures/Figure_6'
plt.savefig(os.path.join(savepath, 
                         'layer_quantification_{0}.pdf'.format(fname)),
            bbox_inches='tight', 
            pad_inches=0.3, format='pdf', transparent = True, dpi=300)

    