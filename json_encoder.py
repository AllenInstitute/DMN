# -*- coding: utf-8 -*-
"""
Created on Wed Jul 05 11:05:26 2017

@author: jenniferwh
"""
import json
import numpy as np
import os

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)
          
#jsondat = os.path.join(path, 'cre_ss_root_ctx_hipp.json')
#with open(jsondat, 'r') as data_file:    
#    data = json.load(data_file)
#with open(os.path.join(outpath, 'filename.json'), 'w') as outfile:
#    json.dump(dat, outfile, sort_keys = False, indent = 4, cls = MyEncoder)
