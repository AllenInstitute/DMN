# -*- coding: utf-8 -*-
"""
Created on Tue Feb 07 17:27:14 2017

get_CAV_grids

@author: jenniferwh
"""

import os
import lims_utilities as lu
from anatomy.cav_grid_cache import CavGridCache
from allensdk.config.manifest import Manifest

basepath = r'E:\CAV_grid'

manifest_file = os.path.join(basepath, 'cav_grid_manifest.json')
Manifest.safe_make_parent_dirs(manifest_file)

cgc = CavGridCache(manifest_file = manifest_file)

image_series_query = '''
    select im.id from image_series im
    join specimens_workflows sw on sw.specimen_id = im.specimen_id
    where sw.workflow_id = 471789262 
    and im.workflow_state like 'passed'
'''
image_series_results = lu.query(image_series_query)
print len(image_series_results)

isr = [iser['id'] for iser in image_series_results]

for im in isr:
    print im
    try:
        cgc.get_cav_grid(im)
    except:
        continue
    
#workflow id = 471789262 for T503.2, 304950893 for T601.4
    