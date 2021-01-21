# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:13:06 2020

@author: jenniferwh
"""
import os
import numpy as np
import nrrd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import h5py
import argparse
import json
from anatomy.anatomy_api import AnatomyApi
import platform

def calculate_injection_centroid(injection_density,
                                 injection_fraction,
                                 resolution=100):
    '''
    Compute the centroid of an injection site.
    
    Parameters
    ----------
    
    injection_density: np.ndarray
        The injection density volume of an experiment

    injection_fraction: np.ndarray
        The injection fraction volume of an experiment

    '''

    # find all voxels with injection_fraction > 0
    injection_voxels = np.nonzero(injection_fraction)
    injection_density_computed = np.multiply(injection_density[injection_voxels],
                                             injection_fraction[injection_voxels]) 
    sum_density = np.sum(injection_density_computed)

    # compute centroid in CCF coordinates
    if sum_density > 0 :
        centroid = np.dot(injection_density_computed,
                          list(zip(*injection_voxels))) / sum_density * resolution
    else:
        centroid = None
    
    return centroid

def get_top_centroid_coord(x,y,z,volume_shape,flipped, path):
    with h5py.File(os.path.join(path, 'surface_paths_10.h5'), 'r') as f:
        path_lookup = f['volume lookup'][()]
        paths = f['paths'][()]
    with h5py.File(os.path.join(path, 'top.h5'), 'r') as f:
        top_lookup = f['view lookup'][()]
        top_size = f.attrs['view size']

    path_indices = path_lookup.flat[top_lookup[:,1]]
    
    top_paths = paths[path_indices,:]
    cortical_proj_blank = np.zeros(top_size)  
    radius = 1
    point_blank = np.zeros(volume_shape)
    point_blank[x-radius:x+radius*2,y-radius:y+radius*2,z-radius:z+radius*2] = 1
    print('first time')
    print(radius)
    if flipped:
        point_blank = np.flip(point_blank,axis=2)
    for streamline in np.arange(top_paths.shape[0]):
        cortical_proj_blank.flat[top_lookup[streamline,0]] = \
                max(point_blank.flat[top_paths[streamline,:]])
    coords = np.round(np.mean(np.where(cortical_proj_blank),axis=1)).astype(int)
    print(coords)
    while coords[0] < 0:
        if radius < 10:
            radius += 1
        else:
            radius *= 10
        point_blank = np.zeros(volume_shape)
        point_blank[x-radius:x+radius*2,y-radius:y+radius*2,z-radius:z+radius*2] = 1
        if flipped:
            point_blank = np.flip(point_blank,axis=2)
        for streamline in np.arange(top_paths.shape[0]):
            cortical_proj_blank.flat[top_lookup[streamline,0]] = \
                max(point_blank.flat[top_paths[streamline,:]])
        coords = np.round(np.mean(np.where(cortical_proj_blank),axis=1)).astype(int)
        print(radius)
        print(coords)
    return coords

def get_flat_centroid_coord(x,y,z,volume_shape,flipped, path):
    with h5py.File(os.path.join(path, 'surface_paths_10.h5'), 'r') as f:
        path_lookup = f['volume lookup'][()]
        paths = f['paths'][()]

    #Flat map lookup table
    flatmap_lookup = h5py.File(os.path.join(path, 'flatmap_lookup.h5'))
    flat_paths = flatmap_lookup['flat_paths'][()]
    flat_map = flatmap_lookup['flat_map'][()]
    fsize = flatmap_lookup.attrs['size']
    flatmap_lookup.close()
    
    path_indeces = path_lookup.flat[flat_paths]
    flat_paths = paths[path_indeces,:]
    flat_blank = np.zeros(fsize,dtype=np.uint8)
    print(flat_blank.shape)
    
    radius = 1
    point_blank = np.zeros(volume_shape,dtype=np.uint8)
    point_blank[x-radius:x+radius*2,y-radius:y+radius*2,z-radius:z+radius*2] = 1
    if flipped:
        point_blank = np.flip(point_blank,axis=2)
    for streamline in np.arange(flat_paths.shape[0]):
        flat_blank.flat[flat_map[streamline]] = max(point_blank.flat[flat_paths[streamline,:]])
    coords = np.round(np.mean(np.where(flat_blank),axis=1)).astype(int)
    print('first flat')
    print(coords)
    while coords[0] < 0:
        if radius < 10:
            radius += 1
        else:
            radius *= 10
        point_blank = np.zeros(volume_shape,dtype=np.uint8)
        point_blank[x-radius:x+radius*2,y-radius:y+radius*2,z-radius:z+radius*2] = 1
        if flipped:
            point_blank = np.flip(point_blank,axis=2)
        for streamline in np.arange(flat_paths.shape[0]):
            flat_blank.flat[flat_map[streamline]] = max(point_blank.flat[flat_paths[streamline,:]])
        coords = np.round(np.mean(np.where(flat_blank),axis=1)).astype(int)
        print(radius)
        print(coords)
    return coords
    
def main():
    if platform.system() == 'Windows':
        path = r'cluster_code\cortical_projections'
    else:
        path = r'cluster_code/cortical_projections'
    # Load reference space
    mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                                 resolution=10)
    avg,meta = mcc.get_template_volume()
    print(avg.shape)
    
    # find centroid in CCF space
    mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json',
                                 resolution=25)
    aapi = AnatomyApi()
    try:
        INJDEN, _ = mcc.get_injection_density(experiment_id=args.id)
        INJFRAC, _ = mcc.get_injection_fraction(experiment_id=args.id)
        DATA_MASK, _ = mcc.get_data_mask(experiment_id=args.id)
    except:
        storage_directory = aapi.get_storage_directory(args.id)
        INJDEN, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                     'injection_density_25.nrrd'))
        INJFRAC, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                        'injection_fraction_25.nrrd'))
        DATA_MASK, _ = nrrd.read(os.path.join(storage_directory, 'grid',
                                           'data_mask_25.nrrd'))
    inj_den = np.multiply(INJDEN, DATA_MASK)
    centroid = calculate_injection_centroid(inj_den, INJFRAC, 25)
    flipped = False
    if centroid[2] > 5700:
        flipped = True
    x,y,z = (centroid/10).astype(int)
    print('get top centroid')
    top_centroid = get_top_centroid_coord(x,y,z,avg.shape,flipped, path)
    print('get flat centroid')
    flat_centroid = get_flat_centroid_coord(x,y,z,avg.shape,flipped, path)
    centroids = {'id': args.id,
                 'top_x': float(top_centroid[0]),
                 'top_y': float(top_centroid[1]),
                 'flat_x': float(flat_centroid[0]),
                 'flat_y': float(flat_centroid[1])}   
    with open(os.path.join(path, 
                           'output', 
                           'cortical_projection_coordinates',
                           '{}.json'.format(args.id)), 'w') as outfile:
        json.dump(centroids, outfile, sort_keys = False, indent = 4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('id', type = int, help='image series id')
    args = parser.parse_args()

    main()
