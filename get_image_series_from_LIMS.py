# -*- coding: utf-8 -*-
"""
Created on Tue Feb 07 17:27:14 2017

get_CAV_grids

@author: jenniferwh
"""

import lims_utilities as lu

#image_series_query = '''
#    select im.id from image_series im
#    join specimens_workflows sw on sw.specimen_id = im.specimen_id
#    where sw.workflow_id = 471789262 
#    and im.workflow_state like 'passed'
#'''

image_series_query = '''
    select im.id from image_series im
    join specimens sp on sp.id = im.specimen_id
    where sp.name like '%Cre%'
    and im.project_id = 24
    and im.workflow_state like 'passed'
'''
    
image_series_results = lu.query(image_series_query)
print len(image_series_results)

isr = [iser['id'] for iser in image_series_results]
print(isr)
    
#workflow id = 471789262 for T503.2, 304950893 for T601.4