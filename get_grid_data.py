# -*- coding: utf-8 -*-
"""
Created on Thu May 23 20:35:57 2019

@author: jenniferwh
"""
import os
import nrrd
import pandas as pd
import json

from anatomy.anatomy_api import AnatomyApi
aapi = AnatomyApi()

isids = pd.read_csv(r"C:\Users\jenniferwh\Desktop\isids.csv")
savepath = r"C:\Users\jenniferwh\Dropbox (Allen Institute)\AD_grid_files"

isids = isids['isids'].values

not_found = []
for isid in isids:
    try:
        gridpath = os.path.join(aapi.get_storage_directory(isid), 'grid')
        print(gridpath)
        injdens, _ = nrrd.read(os.path.join(gridpath, 'injection_density_100.nrrd'))
        injfrac, _ = nrrd.read(os.path.join(gridpath, 'injection_fraction_100.nrrd'))
        projdens, _ = nrrd.read(os.path.join(gridpath, 'projection_density_100.nrrd'))
        outpath = os.path.join(savepath, str(isid))
        if not os.path.isdir(outpath):
            os.mkdir(outpath)
        nrrd.write(os.path.join(outpath, 'injection_density_100.nrrd'), injdens)
        nrrd.write(os.path.join(outpath, 'injection_fraction_100.nrrd'), injfrac)
        nrrd.write(os.path.join(outpath, 'projection_density_100.nrrd'), projdens)
    except:
        not_found += isid
        
jsondat = {'image_series_id': not_found}
with open(os.path.join(savepath, 'not_found.json'), 'w') as outfile:
    json.dump(jsondat, outfile, sort_keys = False, indent = 4)
