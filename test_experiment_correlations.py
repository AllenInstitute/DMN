# tests for anatomy.utilities
from __future__ import division, print_function, absolute_import
import os

import numpy as np
import pytest
import nrrd

from anatomy import utilities


def test_isometric_block_downsample():
    
    image = np.zeros((10, 10))
    
    image[0, 0] = 1
    image[1, 1] = 1
    image[8, 9] = 1
    
    reduced = utilities.isometric_block_downsample(image, 3)
    
    assert( reduced.shape[0] == 4 ) # padding
    assert( reduced[0, 0] == 2 )
    assert( reduced[2, 3] == 1 )
