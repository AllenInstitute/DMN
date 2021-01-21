# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 19:41:58 2017
Nile already generated json files containing the data for intersecting each 
structure with the fMRI overlay map. Here I am summing the overlap of each 
structure with fMRI map values of 0, 1, 2, 3, and 4 and summarizing the data 
for all structures.

@author: jenniferwh
"""

import os
import nrrd
import numpy as np
import pandas as pd

from allensdk.api.queries.ontologies_api import OntologiesApi

#%%
def desc_map(): #Descendant map for easily converting between parents and descendants
    oapi = OntologiesApi()
    sts = oapi.get_structures(1)
    
    dmap = {}
    for st in sts:
        stid = st['id']
        dmap[stid] = map(lambda y: y['id'], filter(lambda x: str(stid) in x['structure_id_path'].split('/'), sts))
    
    return dmap

def find_parents(structures): #Finds parent structure ids 
                                    #for structures at a lower ontology level than summary structures
    dmap = desc_map()
    summary_structures = get_ssids()
    
    parents = []
    no_parent = [] #keep track of structures with a higher ontology level
    
    for st in structures:
        if st in summary_structures:
            parents.append(st)
        else:
            keys = []
            for key, val in dmap.items():
                if st in val:
                    keys.append(key)
            parent = [ k for k in keys if k in summary_structures ]
            if parent == []:
                parent = st
                no_parent.append(st)
            else:
                parent = parent[0]
            parents.append(parent)
    return parents, no_parent

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

def get_ssids():
    oapi = OntologiesApi()
    return [st['id'] for st in oapi.get_structures(structure_set_names = "'Brain – Summary Structures'")]  

def get_ssnames():
    oapi = OntologiesApi()
    return [[st['id'], st['name']] for st in oapi.get_structures(structure_set_names = "'Brain – Summary Structures'")]  
    
#%%
basepath = r'2019 DMN\Data files'
intersect_path = r'DMN\structure_intersections'

#%%
## Iterate through structure intersect files and sum the overlapping region to find the total size of overlap

structures = []
fmrivals = []
sumpix = []

for fname in os.listdir(intersect_path):
    structures.append(fname[10:-12])
    fmrivals.append(fname[-6])
    dat, meta = nrrd.read(os.path.join(intersect_path, fname))
    sumpix.append(np.sum(np.nonzero(dat)))

stids = [int(float(structure)) for structure in structures]
fmrivals = [int(float(val)) for val in fmrivals]
sumpix = [int(float(pix)) for pix in sumpix]

ss_ids, leftovers = find_parents(stids)

str_names = find_names(ss_ids, 'name')

#%%
## Check the identity of structures that are part of the fMRI mask but were not found in the set of summary structures. 
## These structures will still be quantified for their overlap wth the fMRI mask.

unique_leftovers = set(leftovers)
print(find_names(unique_leftovers, 'name'))

#%%
## Create a pandas dataframe to hold all the data
## Reshape arrays so that there will be one row for each region and columns for each fMRI value
ssids_res = np.reshape(ss_ids, (int(len(ss_ids)/3),3))
fmrivals_res = np.reshape(fmrivals, (int(len(fmrivals)/3),3))
sumpix_res = np.reshape(sumpix, (int(len(sumpix)/3),3))
names_res = np.reshape(str_names, (int(len(str_names)/3),3))
acronyms = find_names(ssids_res[:,0], 'acronym')
df = pd.DataFrame([ssids_res[:,0],
                   ssids_res[:,0],
                   names_res[:,0],
                   acronyms,
                   fmrivals_res[0], 
                   fmrivals_res[1],
                   fmrivals_res[2]],
                  index = ['structure_id', 'stid', 'name', 'acronym', '0', '1', '2'])
df = df.transpose()
df['0'] = sumpix_res[:,0]
df['1'] = sumpix_res[:,1]
df['2'] = sumpix_res[:,2]
df = df.groupby(['structure_id', 'name','acronym']).sum().reset_index()
#%%
##Count total voxels and calculate percent of each region in DMN masks
df['total_voxels'] = sum([df['0'], df['1'], df['2']])
df['percent1'] = np.round(sum([df['1']/df['total_voxels']])*100, 2)
df['percent2'] = np.round(sum([df['2']/df['total_voxels']])*100, 2)
df['percent12'] = np.round((sum([df['1'], df['2']])/df['total_voxels'])*100, 2)

#%%
#SAVE
df.to_csv(os.path.join(basepath, 'summary_structures_percent_overlap.csv'),
          index = False)
