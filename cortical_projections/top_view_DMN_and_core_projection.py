#!/shared/utils.x86_64/python-2.7/bin/python
import os
import nrrd
import numpy as np
import h5py
from PIL import Image
import argparse
import scipy.ndimage
import matplotlib.pyplot as plt

MODEL_DIRECTORY = '/projects/0378/vol1/informatics/model/P56'

def get_paths_lut( view_file ):

    # open view index file
    vi = h5py.File( view_file, 'r' )
    
    # initialize
    lut = vi['view lookup'][:]
    paths = vi['paths'][:]
    
    vi.close()
    
    return {'lut': lut, 'paths': paths}


def create_top_view_projection( image, top_view, output_default_value = 0 ):

    # initialize
    output = np.zeros( top_view['lut'].shape + (2,), image.dtype )
    output.fill( output_default_value )
    output[:,:,1] = 0
    
    def max_along_path( pid ) :
        path = top_view['paths'][pid,:]
        arr = image.flat[ path ]
        idmax = np.argmax( arr )
        vmax = arr[ idmax ]
        indmax = path[ idmax ]

        return vmax, indmax
    
    ind = np.where( top_view['lut'] > -1 )
    output[ind] = map( max_along_path, top_view['lut'][ind] )
    
    return output[:,:,0], output[:,:,1].astype(np.int32)


def compute_depth_coordinates( depth, template_shape, top_view, default_value=-1 ):
    xyz = np.unravel_index(depth[:], template_shape)

    bgind = np.where( top_view['lut'] < 0 )

    x = xyz[0].reshape(depth.shape).copy()
    y = xyz[1].reshape(depth.shape).copy()
    z = xyz[2].reshape(depth.shape).copy()
    
    x[bgind] = default_value
    y[bgind] = default_value
    z[bgind] = default_value

    coords = np.zeros( (3,) + depth.shape, dtype = np.int32 )
    coords[0,:,:] = x
    coords[1,:,:] = y
    coords[2,:,:] = z

    return coords


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_directory', default=MODEL_DIRECTORY)
    parser.add_argument('--resolution', type=int, default=100)
    parser.add_argument('output_file')
    args = parser.parse_args()


    # hardcoding this for now
    args.resolution = 10

    template_file = '/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/cluster_code/cortical_projections/dmn_mask_and_core.nrrd'
    
    view_file = os.path.join( args.model_directory, 'corticalCoordinates', 'top_view_paths_%d.h5' % args.resolution)

    top_view = get_paths_lut( view_file )

    template_100, meta = nrrd.read( template_file )
    dmn_100 = np.zeros(template_100.shape)
    dmn_100[np.where(template_100 > 0)] = 1
    core_100 = np.zeros(template_100.shape)
    core_100[np.where(template_100 == 2)] = 1
    dmn = scipy.ndimage.zoom(dmn_100, 10, order=4)
    dmn[(0,0,0)] = 0
    core = scipy.ndimage.zoom(core_100, 10, order=4)
    core[(0,0,0)] = 0

    dmn_output, depth = create_top_view_projection( dmn, top_view )
    ds_dmn_output = scipy.ndimage.zoom(dmn_output, 0.4, order=4)
    core_output, depth = create_top_view_projection( core, top_view )
    ds_core_output = scipy.ndimage.zoom(core_output, 0.4, order=4)
    
    nrrd.write( os.path.join(args.output_file, 'DMN_mask_cortical_projection_10.nrrd'), dmn_output )
    nrrd.write( os.path.join(args.output_file, 'DMN_mask_cortical_projection_25.nrrd'), ds_dmn_output )
    nrrd.write( os.path.join(args.output_file, 'core_mask_cortical_projection_10.nrrd'), core_output )
    nrrd.write( os.path.join(args.output_file, 'core_mask_cortical_projection_25.nrrd'), ds_core_output )
    im = Image.fromarray(np.uint8(plt.cm.gray(ds_dmn_output))*255)
    im.save(os.path.join(args.output_file, 'DMN_cortical_projection_25.png'))
    im = Image.fromarray(np.uint8(plt.cm.gray(ds_core_output))*255)
    im.save(os.path.join(args.output_file, 'core_cortical_projection_25.png'))
    im = Image.fromarray(np.uint8(plt.cm.gray(dmn_output))*255)
    im.save(os.path.join(args.output_file, 'DMN_cortical_projection_10.png'))
    im = Image.fromarray(np.uint8(plt.cm.gray(core_output))*255)
    im.save(os.path.join(args.output_file, 'core_cortical_projection_10.png'))


if __name__ == '__main__': main()


