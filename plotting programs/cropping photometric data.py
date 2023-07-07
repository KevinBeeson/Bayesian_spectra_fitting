#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:08:22 2023

@author: kevin
"""
from astropy.io.votable import parse
import numpy as np
from astropy.io.votable import from_table,writeto

votable = parse("Ruprecht_147_photometric_cross.xml")
photometric_data=votable.get_first_table().to_table(use_names_over_ids=True)
data=np.loadtxt("minimal_run.txt",dtype=np.int64)
mask=[x['sobject_id'] in data for x in photometric_data ]
photometric_data_small=photometric_data[mask]
votable=from_table(photometric_data_small)
writeto(votable,"minimum_photometric_cross.xml")
