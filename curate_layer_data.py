# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 19:37:48 2020

@author: jenniferwh
"""

import os
import pandas as pd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from anatomy.anatomy_api import AnatomyApi

import platform
if platform.system() == 'Windows':
    basepath = r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Mesoscale Connectome Papers in Progress\2019 DMN'
elif platform.system() == 'Darwin':
    basepath = r'/Users/jenniferwh/Dropbox (Allen Institute)/Mesoscale Connectome Papers in Progress/2019 DMN'

datapath = os.path.join(basepath, 'data_files')
savepath = os.path.join(basepath, '_new_figures', 'Figure_3')

dat = pd.read_csv(os.path.join(savepath, 'Ext Data Table 3_Corticocortical NPV matrices.csv'))
dat.rename(columns = {'Exp image series ID': 'image_series_id', 
                      'Rbp4 anchor image series ID': 'Rbp4_anchor_id',
                     'distance from anchor': 'distance', 'Consensus Rbp4 anchor Source': 'Rbp4_source'}, 
           inplace = True)
dmn_dat = pd.read_csv(os.path.join(datapath, 
                                   'wt_cre_matched_layer_injections_DMN_and_core_projections_coefficients.csv'))
dmn_dat.rename(columns = {'id': 'image_series_id'}, inplace = True)

dat = dat.merge(dmn_dat[['image_series_id', 'injection dmn fraction',
                   'projection dmn fraction', 
                   'distance coefficient', 
                   'DMN coefficient',
                        'injection core fraction',
                        'projection core fraction',
                        'core distance coefficient',
                        'DMN core coefficient']], 
          on = 'image_series_id', how = 'left')
dat.loc[dat['Mouse Line'] == 'C57BL/6J / Emx1', 'layer'] = 'all'
dat.loc[dat['Mouse Line'].isin(['Cux2-IRES-Cre', 'Sepw1-Cre_NP39']), 'layer'] = '2/3'
dat.loc[dat['Mouse Line'].isin(['Nr5a1-Cre', 'Scnn1a-Tg3-Cre', 'Rorb-IRES2-Cre']), 'layer'] = '4'
dat.loc[dat['Mouse Line'].isin(['Tlx3-Cre_PL56']), 'layer'] = '5 IT' #separate
dat.loc[dat['Mouse Line'].isin(['Rbp4-Cre_KL100']), 'layer'] = '5 IT PT' #separate
dat.loc[dat['Mouse Line'].isin(['Chrna2-Cre_OE25', 'Efr3a-Cre_NO108', 'Sim1-Cre_KJ18',
                               'A930038C07Rik-Tg1-Cre']), 'layer'] = '5 PT'
dat.loc[dat['Mouse Line'].isin(['Ntsr1-Cre_GN220', 'Syt6-Cre_KI148']), 'layer'] = '6'
# make each source have only in or only out injections
drop_expts = [477037203, #single MOs injection outside the DMN
              141602484, #another MOs injection outside the DMN (barely: 46%)
             300929973, #one VISp injection inside the DMN
             272821309] #another VISp injection inside the DMN
drop_expts += list(dat[(dat['Rbp4_source'] == 'VISp-2') & 
                      (dat['injection dmn fraction'] > 0.5)]['image_series_id'].values)
drop_expts += list(dat[(dat['Rbp4_source'] == 'SSp-bfd') & 
                      (dat['injection dmn fraction'] < 0.5)]['image_series_id'].values)
drop_expts += list(dat[dat['Rbp4_source'] == 'VISpm']['image_series_id'].values)
dat = dat[~dat['image_series_id'].isin(drop_expts)]
drop_5PT = [278258073, 267750528, 156394513, 176433237, 287807030, 287995889]
dat = dat[~dat['image_series_id'].isin(drop_5PT)]
drop_L6 = dat[(dat['layer'] == '6') &
             (dat['Manual PN class '] != 'CT')]['image_series_id'].values
dat = dat[~dat['image_series_id'].isin(drop_L6)]
# Only keep PT experiments for L5 PT
drop_L5 = dat[(dat['layer'] == '5 PT') &
             (~dat['Manual PN class '].isin(['PT', 'local']))]['image_series_id'].values
dat = dat[~dat['image_series_id'].isin(drop_L5)]
# One Rbp4 experiment was classified as PT
drop_L5 = dat[(dat['layer'] == '5 IT PT') &
             (dat['Manual PN class '] == 'PT')]['image_series_id'].values
dat = dat[~dat['image_series_id'].isin(drop_L5)]
print(len(dat))
dat.to_csv(os.path.join(savepath, 'curated_cre_layer_data.csv'), index = False)