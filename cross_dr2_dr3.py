#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 15:35:34 2022

@author: kevin
"""
import csv
from astropy.io.votable import parse
from astropy.io import fits
import numpy as np

votable=parse('all_gaia.xml')
gaia=votable.get_first_table().to_table(use_names_over_ids=True)
votable = parse("1642685188743O-result.vot")
dr2_dr3=votable.get_first_table().to_table(use_names_over_ids=True)
dr2_dr3=np.vstack([dr2_dr3['dr3_source_id'],dr2_dr3['dr2_source_id']])
# a=np.array(gaia['source_id'],dtype=str)
# np.savetxt('cross_match',a,fmt='%s')
# gaia_cross = np.array(list(csv.reader(open('gaia_cross/skymapper2BestNeighbour0001.csv', 'rt'), delimiter=',')))

# hermes=fits.open('dr6.0.fits')