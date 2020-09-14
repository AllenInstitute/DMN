# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 20:15:29 2020

@author: jenniferwh
"""

import os
import json
import pandas as pd

savepath = r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\cluster_code\cortical_projections\output'
readpath = r'\\allen\programs\celltypes\workgroups\mousecelltypes\T503_Connectivity_in_Alzheimer_Mice\Jennifer\cluster_code\cortical_projections\output\cortical_projection_coordinates'

df = pd.DataFrame(columns = {'id', 'top_x', 'top_y', 'flat_x', 'flat_y'})
fnames = []
for filename in os.listdir(readpath):
    print(filename)
    jsondat = os.path.join(readpath, filename)
    with open(jsondat, 'r') as data_file:    
        dat = json.load(data_file)
    df = df.append(dat, ignore_index = True)    

df.to_csv(os.path.join(savepath, 'cortical_flatmap_coordinates.csv'),
          index = False)