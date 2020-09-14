# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:13:06 2020

@author: jenniferwh
"""
import os
import pandas as pd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import h5py
import argparse
import json
    
def get_top_centroid_coord(x,y,z,volume_shape,radius):
    path = r'/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/cluster_code/cortical_projections'
    with h5py.File(os.path.join(path, 'surface_paths_10.h5'), 'r') as f:
        path_lookup = f['volume lookup'][()]
        paths = f['paths'][()]
    with h5py.File(os.path.join(path, 'top.h5'), 'r') as f:
        top_lookup = f['view lookup'][()]
        top_size = f.attrs['view size']
    # This gives us the path indices that point to the side
    path_indices = path_lookup.flat[top_lookup[:,1]]
    # Now we can subsample the paths
    top_paths = paths[path_indices,:]
    point_blank = np.zeros(volume_shape)
    flat_blank = np.zeros(top_size)    
    point_blank[x-radius:x+radius,y-radius:y+radius,z-radius:z+radius] = 1
    
    for streamline in np.arange(top_paths.shape[0]):
        flat_blank.flat[top_lookup[streamline,0]] = \
                max(point_blank.flat[top_paths[streamline,:]])
    
    return np.round(np.mean(np.where(flat_blank),axis=1)).astype(int)

def main():
    path = r'/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/cluster_code/cortical_projections'
    mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                                 resolution=10)
    avg,meta = mcc.get_template_volume()
    structure_tree = mcc.get_structure_tree()
    iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
    ctx_exps = pd.DataFrame(mcc.get_experiments(cre=False, 
                                           injection_structure_ids=[iso['id']]))
    cre_exps = pd.DataFrame(mcc.get_experiments(cre=['Emx1-IRES-Cre','Rbp4-Cre_KL100'],
                                          injection_structure_ids = [iso['id']]))
    ctx_exps = pd.concat([ctx_exps, cre_exps])
    fail_expts = [114008926, 120280939, 180073473, 180403712, 180601025, 183174303, 183329222,
                  249396394, 296047806, 299446445, 301060890, 303784745, 480069939, 482578964, 
                  506947040, 514333422, 525796603, 545428296, 559878074, 638314843, 182888003,
                 304585910, 183171679, 272930013, 523718075, 517072832, 148964212, 304762965,
                 566992832, 272930013, 304762965, 266250904, 114399224, 286483411, 286417464,
                 593277684, 546103149, 642809043, 286483411, 304564721] #VISp outlier excluded
    ctx_exps = ctx_exps[~ctx_exps['id'].isin(fail_expts)]
    subset = ctx_exps[ctx_exps['id'] == args.image_series_id]
    # Put all injections in left hemisphere but keep track of where they started
    flipped = False
    if subset['injection_z'].values[0] > 5700:
        flipped = True
        subset.loc[subset['id'] == args.image_series_id, 'injection_z'] = subset['injection_z'] - 5700
    x,y,z = (subset[['injection_x','injection_y','injection_z']]/10).astype(int).values[0]
    top_centroid = get_top_centroid_coord(x,y,z,avg.shape,1)
    if top_centroid[0] < 0:
        j = 1
        while top_centroid[0] < 0:
            j += 1
            if j > 10:
                top_centroid = get_top_centroid_coord(x,y,z,avg.shape,j*10)
            else:
                top_centroid = get_top_centroid_coord(x,y,z,avg.shape,j)
    if flipped:
        top_centroid[0] = 570+(570-top_centroid[0])
    centroids = {'id': args.image_series_id,
                 'top_x': float(top_centroid[0]),
                 'top_y': float(top_centroid[1])}   
    with open(os.path.join(path, 'output', '{}.json'.format(args.image_series_id)), 'w') as outfile:
        json.dump(centroids, outfile, sort_keys = False, indent = 4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('image_series_id', type = int, help='image series id')
    args = parser.parse_args()

    main()