# -*- coding: utf-8 -*-
"""
Created on Tue Feb 07 17:27:14 2017

Get image series ids only

@author: jenniferwh
"""

import lims_utilities as lu
from itertools import chain

def image_series_query(experiment_type):
    if experiment_type == 'target_defined':
        workflows = [471789262, 304950893]
    image_series_ids = []
    for workflow in workflows:
        image_series_query = '''
            select im.id, sp.name from image_series im
            join specimens sp on sp.id = im.specimen_id
            join specimens_workflows sw on sw.specimen_id = im.specimen_id
            where sw.workflow_id = %s
            and im.workflow_state like 'passed'
        '''%workflow

        image_series_results = lu.query(image_series_query)
        print len(image_series_results)

        isr = [iser['id'] for iser in image_series_results]
        image_series_ids.append(isr)

    image_series_list = list(chain.from_iterable(image_series_ids))
    return image_series_list

# image_series_query = '''
#     select im.id from image_series im
#     join specimens sp on sp.id = im.specimen_id
#     where sp.name like '%Cre%'
#     and im.project_id = 24
#     and im.workflow_state like 'passed'
# '''
#workflow id = 471789262 for T503.2, 304950893 for T601.4