# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 17:50:45 2020

@author: jenniferwh
"""
import pytest

def test_r_hem_experiment():
    import get_cortical_projection_coords

    assert get_cortical_projection_coords.main(307297141) == [911 956]
