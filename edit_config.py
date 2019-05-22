# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 13:40:40 2017

@author: jenniferwh
"""

import os
import json
from json_encoder import MyEncoder
from anatomy.anatomy_api import AnatomyApi

path = r'E:\fMRI_unionize\inputs\configs'

jsondat = os.path.join(path, 'aims_2_and_4_ss_root_ctx_hipp.json')
with open(jsondat, 'r') as data_file:    
    data = json.load(data_file)

aapi = AnatomyApi()
ims = aapi.get_image_series_by_workflow([471789262, 304950893])
data['image_series_ids'] = ims
data['plaque_dataset'] = False

with open(os.path.join(path, 'aims_2_and_4_ss_root_ctx_hipp.json'), 'w') as outfile:
    json.dump(data, outfile, sort_keys = False, indent = 4, cls = MyEncoder)