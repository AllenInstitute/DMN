#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 19:29:04 2020

@author: jenniferwh
"""

import cortical_flatmap_projection
import pandas as pd

dat = pd.read_csv(r'/Users/jenniferwh/Dropbox (Personal)/2019 DMN/python_code/DMN/data_files/cortical_flatmap_coordinates.csv')
isids = dat['id'].values
for isid in isids:
    print(isid)
    cortical_flatmap_projection.main(isid)