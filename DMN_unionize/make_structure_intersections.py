#!/shared/utils.x86_64/python-2.7/bin/python

#General imports
from __future__ import division, print_function, absolute_import
import os
import numpy as np
import pandas as pd
from fmri_intersection import FmriIntersect
import platform
import nrrd
        
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache()
st = mcc.get_structure_tree()

def make_structure_intersections(maskpath, force_update = False, structures = None):
    ''' Find overlap between CCF structures and a supplied mask
    Parameters
    ----------
    maskpath : str
        path to a ccf-shaped mask as an nrrd file
    force_update : bool, optional
        if True -> CCF structure mask will be computed 
        regardless of whether it can be loaded. Defaults 
        to False
        
    Returns
    -------
    
    '''
    basepath = os.getcwd()
    fmri_values = [0, 1, 2]
    
    inter = FmriIntersect(
            top_level=basepath, 
            fmri_path=maskpath, 
            )
    inter.make_fmri_masks()
    
    if structures == None:
        oapi = OntologiesApi()
        structs = pd.DataFrame(oapi.get_structures(structure_set_names="'Mouse Connectivity - Summary'"))
        new_structs = st.get_structures_by_acronym(['ProS', 'HATA', 'APr'])
        structure_ids = list(structs['id'].values)
        structure_ids.extend([struct['id'] for struct in new_structs])
        structure_ids.remove(934)
    else:
        structure_ids = structures
    inter.many_structure_masks(structure_ids, force_update=True)
    inter.many_intersection_masks(structure_ids, fmri_values)

def find_names(structures, kind = 'name'): #Finds name or acronym for list of structure ids.
                                            # Queries entire set of mouse structures, not just summary structures
    oapi = OntologiesApi()
    all_structures = oapi.get_structures(1)
    
    strnames = []
    for st in structures:
        for struct in all_structures:
            if struct['id'] == st:
                if kind == 'name':
                    stname = struct['safe_name']
                elif kind == 'acronym':
                    stname = struct['acronym']
        strnames.append(stname)
    return strnames

if platform.system() == 'Darwin':      
    path = r'dmn_mask_and_core.nrrd'
elif platform.system() == 'Windows':
    path = r'dmn_mask_and_core.nrrd'
make_structure_intersections(path)
